import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go


st.set_page_config(
    page_title="CysSelect",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

APP_NAME = "CysSelect"
APP_TAGLINE = (
    "Chemoproteomic hit calling, enantiomer comparison, prioritization, "
    "simulation, and next-library design for covalent fragment discovery."
)

DEFAULT_CONFIG = {
    "liganded_log2cr_threshold": 1.5,
    "enantioselective_abs_delta_auc_threshold": 15.0,
    "sweet_spot_min_sites": 4,
    "sweet_spot_max_sites": 64,
    "sweet_spot_min_fraction": 0.25,
    "pvalue_cutoff": 0.05,
    "n_sim_pairs": 73,
    "n_sim_cysteines": 700,
    "sim_noise_sd": 0.20,
    "sim_missing_rate": 0.08,
    "seed": 42,
}

REQUIRED_COLUMNS_REAL = [
    "pair_id",
    "enantiomer",
    "cysteine_id",
    "concentration_uM",
    "replicate",
    "intensity_dmso",
    "intensity_condition",
]

OPTIONAL_COLUMNS_REAL = [
    "p_value",
    "protein_id",
    "gene_name",
    "scaffold_cluster",
]


def init_state() -> None:
    defaults = {
        "config": DEFAULT_CONFIG.copy(),
        "raw_df": None,
        "processed_hits_df": None,
        "pair_summary_df": None,
        "sim_raw_df": None,
        "sim_hits_df": None,
        "sim_pair_summary_df": None,
        "virtual_library_df": None,
        "scored_library_df": None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def reset_session() -> None:
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    init_state()


def demo_real_dataset() -> pd.DataFrame:
    rng = np.random.default_rng(42)
    rows = []
    pairs = [f"pair_{i:03d}" for i in range(1, 16)]
    cysteines = [f"CYS_{i:04d}" for i in range(1, 61)]
    concentrations = [50.0, 200.0]

    for pair in pairs:
        scaffold = int(rng.integers(1, 6))
        for enantiomer in ["R", "S"]:
            stereo_bias = rng.normal(0.0, 0.3)
            for cyst in rng.choice(cysteines, size=18, replace=False):
                base_strength = rng.normal(0.0, 1.0) + (0.4 if enantiomer == "R" else -0.4) * stereo_bias
                for conc in concentrations:
                    for rep in [1, 2, 3]:
                        intensity_dmso = float(np.exp(rng.normal(8.0, 0.3)))
                        effect = 1 / (1 + np.exp(-(base_strength + 0.9 * np.log10(conc / 50.0))))
                        intensity_condition = intensity_dmso * max(0.02, 1 - 0.75 * effect)
                        intensity_condition *= np.exp(rng.normal(0.0, 0.12))

                        rows.append(
                            {
                                "pair_id": pair,
                                "enantiomer": enantiomer,
                                "cysteine_id": cyst,
                                "protein_id": f"P{rng.integers(10000, 99999)}",
                                "gene_name": f"GENE_{rng.integers(1, 300)}",
                                "concentration_uM": conc,
                                "replicate": rep,
                                "intensity_dmso": intensity_dmso,
                                "intensity_condition": intensity_condition,
                                "p_value": float(np.clip(rng.uniform(0.0001, 0.2), 0.0001, 1.0)),
                                "scaffold_cluster": scaffold,
                            }
                        )
    return pd.DataFrame(rows)


def safe_download_button(df: pd.DataFrame, label: str, filename: str, key: str) -> None:
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=label,
        data=csv_bytes,
        file_name=filename,
        mime="text/csv",
        key=key,
    )


def validate_real_data(df: pd.DataFrame) -> dict:
    warnings = []
    errors = []

    for col in REQUIRED_COLUMNS_REAL:
        if col not in df.columns:
            errors.append(f"Missing required column: {col}")

    if errors:
        return {
            "ready": False,
            "warnings": warnings,
            "errors": errors,
            "n_rows": len(df),
            "n_pairs": 0,
            "n_cysteines": 0,
            "missing_fraction": float(df.isna().mean().mean()) if len(df) else 0.0,
        }

    if not set(df["enantiomer"].dropna().astype(str).unique()).issubset({"R", "S"}):
        warnings.append("Enantiomer column contains values other than R and S.")

    pair_counts = df.groupby("pair_id")["enantiomer"].nunique()
    bad_pairs = pair_counts[pair_counts < 2]
    if len(bad_pairs) > 0:
        warnings.append(f"{len(bad_pairs)} pair(s) do not have both R and S entries.")

    if df["concentration_uM"].nunique() < 2:
        warnings.append("Only one concentration detected. AUC and ΔAUC are more informative with two concentrations.")

    return {
        "ready": len(errors) == 0,
        "warnings": warnings,
        "errors": errors,
        "n_rows": len(df),
        "n_pairs": int(df["pair_id"].nunique()),
        "n_cysteines": int(df["cysteine_id"].nunique()),
        "missing_fraction": float(df.isna().mean().mean()),
    }


def compute_cr(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["cr"] = out["intensity_dmso"] / out["intensity_condition"].replace(0, np.nan)
    out["log2_cr"] = np.log2(out["cr"].clip(lower=1e-9))
    out["percent_comp"] = (1.0 - 1.0 / out["cr"].clip(lower=1e-9)) * 100.0
    out["percent_pep"] = (1.0 / out["cr"].clip(lower=1e-9)) * 100.0
    if "p_value" in out.columns:
        out["neglog10_p"] = -np.log10(out["p_value"].clip(lower=1e-300))
    else:
        out["neglog10_p"] = np.nan
    return out


def summarize_hits(
    df: pd.DataFrame,
    liganded_log2cr_threshold: float,
    enantioselective_abs_delta_auc_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = compute_cr(df)

    agg_dict = {
        "cr": ("cr", "mean"),
        "log2_cr": ("log2_cr", "mean"),
        "percent_comp": ("percent_comp", "mean"),
        "percent_pep": ("percent_pep", "mean"),
    }
    if "p_value" in work.columns:
        agg_dict["p_value"] = ("p_value", "mean")
    if "protein_id" in work.columns:
        agg_dict["protein_id"] = ("protein_id", "first")
    if "gene_name" in work.columns:
        agg_dict["gene_name"] = ("gene_name", "first")
    if "scaffold_cluster" in work.columns:
        agg_dict["scaffold_cluster"] = ("scaffold_cluster", "first")

    grp = (
        work.groupby(["pair_id", "enantiomer", "cysteine_id", "concentration_uM"], as_index=False)
        .agg(**agg_dict)
    )

    concentrations = sorted(grp["concentration_uM"].dropna().unique().tolist())
    if len(concentrations) < 2:
        concentrations = [50.0, 200.0]
    c_low = float(concentrations[0])
    c_high = float(concentrations[-1])

    delta_log = np.log10(c_high * 1e-6) - np.log10(c_low * 1e-6)
    if delta_log == 0:
        delta_log = 1.0

    pivot = grp.pivot_table(
        index=["pair_id", "cysteine_id"],
        columns=["enantiomer", "concentration_uM"],
        values=["cr", "log2_cr", "percent_comp", "percent_pep"],
        aggfunc="first",
    )
    pivot.columns = [f"{a}_{b}_{int(c)}" for a, b, c in pivot.columns]
    pivot = pivot.reset_index()

    meta_cols = ["pair_id", "cysteine_id"]
    for c in ["protein_id", "gene_name", "scaffold_cluster"]:
        if c in grp.columns:
            meta_cols.append(c)

    meta = grp[meta_cols].drop_duplicates(subset=["pair_id", "cysteine_id"])
    out = pivot.merge(meta, on=["pair_id", "cysteine_id"], how="left")

    def col_or_nan(name: str) -> pd.Series:
        return out[name] if name in out.columns else pd.Series(np.nan, index=out.index)

    out["auc_r"] = ((col_or_nan(f"percent_comp_R_{int(c_high)}") + col_or_nan(f"percent_comp_R_{int(c_low)}")) / 2.0) * delta_log
    out["auc_s"] = ((col_or_nan(f"percent_comp_S_{int(c_high)}") + col_or_nan(f"percent_comp_S_{int(c_low)}")) / 2.0) * delta_log
    out["delta_auc"] = out["auc_r"] - out["auc_s"]

    cr_candidates = [c for c in out.columns if c.startswith("log2_cr_")]
    out["max_log2_cr"] = out[cr_candidates].max(axis=1)
    out["liganded"] = out["max_log2_cr"] >= liganded_log2cr_threshold

    out["selectivity_class"] = np.where(
        ~out["liganded"],
        "not_liganded",
        np.where(
            out["delta_auc"] >= enantioselective_abs_delta_auc_threshold,
            "R_selective",
            np.where(
                out["delta_auc"] <= -enantioselective_abs_delta_auc_threshold,
                "S_selective",
                "non_enantioselective",
            ),
        ),
    )

    pair_agg = {
        "promiscuity": ("liganded", "sum"),
        "n_total_sites": ("cysteine_id", "nunique"),
        "n_r_selective": ("selectivity_class", lambda s: int((s == "R_selective").sum())),
        "n_s_selective": ("selectivity_class", lambda s: int((s == "S_selective").sum())),
    }
    if "scaffold_cluster" in out.columns:
        pair_agg["scaffold_cluster"] = ("scaffold_cluster", "first")

    pair_summary = out.groupby("pair_id", as_index=False).agg(**pair_agg)
    pair_summary["n_enantioselective"] = pair_summary["n_r_selective"] + pair_summary["n_s_selective"]
    pair_summary["stereoselectivity_fraction"] = np.where(
        pair_summary["promiscuity"] > 0,
        pair_summary["n_enantioselective"] / pair_summary["promiscuity"],
        0.0,
    )

    cfg = st.session_state.config
    pair_summary["sweet_spot"] = (
        (pair_summary["promiscuity"] >= cfg["sweet_spot_min_sites"])
        & (pair_summary["promiscuity"] <= cfg["sweet_spot_max_sites"])
        & (pair_summary["stereoselectivity_fraction"] >= cfg["sweet_spot_min_fraction"])
    )

    return out, pair_summary


def simulate_dataset(
    n_pairs: int,
    n_cysteines: int,
    noise_sd: float,
    missing_rate: float,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    cysteine_ids = [f"SIM_CYS_{i:05d}" for i in range(1, n_cysteines + 1)]
    pair_ids = [f"sim_pair_{i:03d}" for i in range(1, n_pairs + 1)]

    for pair in pair_ids:
        scaffold_cluster = int(rng.integers(1, 20))
        clogp = rng.normal(2.2, 0.9)
        tpsa = max(5.0, rng.normal(55.0, 18.0))
        reactivity = rng.normal(19.5, 1.5)
        stereo_bias = rng.normal(0.0, 0.35)

        selected_cysteines = rng.choice(cysteine_ids, size=int(rng.integers(12, 40)), replace=False)

        for enantiomer in ["R", "S"]:
            handedness = 1.0 if enantiomer == "R" else -1.0

            for cyst in selected_cysteines:
                local_polarity = rng.normal(0.0, 1.0)
                nucleophilicity = rng.normal(0.0, 1.0)
                stereo_sensitivity = rng.normal(0.0, 1.0)

                latent = (
                    0.7 * clogp
                    - 0.35 * reactivity
                    + 0.55 * (tpsa / 100.0) * local_polarity
                    + 0.75 * handedness * stereo_bias * stereo_sensitivity
                    + 0.65 * nucleophilicity
                )

                for concentration in [50.0, 200.0]:
                    for replicate in [1, 2, 3]:
                        dmso = float(np.exp(rng.normal(8.0, 0.25)))
                        signal = latent + 1.2 * np.log10(concentration / 50.0)
                        pct_comp = 100.0 / (1.0 + np.exp(-(signal - 0.5)))
                        cr = 1.0 / max(1e-6, 1 - pct_comp / 100.0)
                        obs_cr = max(1.0, cr * np.exp(rng.normal(0, noise_sd)))

                        intensity_condition = dmso / obs_cr
                        if rng.random() < missing_rate:
                            intensity_condition = np.nan

                        rows.append(
                            {
                                "pair_id": pair,
                                "enantiomer": enantiomer,
                                "cysteine_id": cyst,
                                "protein_id": f"SIM_P{rng.integers(10000, 99999)}",
                                "gene_name": f"SIMGENE_{rng.integers(1, 999)}",
                                "concentration_uM": concentration,
                                "replicate": replicate,
                                "intensity_dmso": dmso,
                                "intensity_condition": intensity_condition,
                                "p_value": float(np.clip(rng.uniform(0.0001, 0.2), 0.0001, 1.0)),
                                "scaffold_cluster": scaffold_cluster,
                                "clogp": clogp,
                                "tpsa": tpsa,
                                "reactivity_proxy": reactivity,
                            }
                        )

    sim_df = pd.DataFrame(rows).dropna(subset=["intensity_condition"]).reset_index(drop=True)
    return sim_df


def generate_virtual_library(n_candidates: int = 500, seed: int = 42) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_candidates):
        rows.append(
            {
                "candidate_id": f"cand_{i:04d}",
                "scaffold_cluster": int(rng.integers(1, 25)),
                "clogp": rng.normal(2.2, 0.9),
                "tpsa": max(5.0, rng.normal(55.0, 18.0)),
                "reactivity_proxy": rng.normal(19.5, 1.5),
                "shape_3d": np.clip(rng.normal(0.5, 0.15), 0.0, 1.0),
            }
        )
    return pd.DataFrame(rows)


def score_virtual_library(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rng = np.random.default_rng(42)

    out["pred_promiscuity"] = np.clip(
        18 + 9 * (out["clogp"] - 2.0) - 2.2 * (out["reactivity_proxy"] - 19.5) + rng.normal(0, 3, len(out)),
        0,
        None,
    )

    out["pred_stereo_fraction"] = np.clip(
        0.35 + 0.004 * (out["tpsa"] - 55.0) - 0.015 * (out["clogp"] - 2.2) + rng.normal(0, 0.06, len(out)),
        0.0,
        1.0,
    )

    def smooth_range_score(x, low=4, high=64, width=6):
        left = 1 / (1 + np.exp(-(x - low) / width))
        right = 1 / (1 + np.exp((x - high) / width))
        return left * right

    def smooth_threshold_score(x, threshold=0.25, width=0.04):
        return 1 / (1 + np.exp(-(x - threshold) / width))

    out["sweet_score_soft"] = (
        0.5 * smooth_range_score(out["pred_promiscuity"], 4, 64, 6)
        + 0.5 * smooth_threshold_score(out["pred_stereo_fraction"], 0.25, 0.04)
    )

    cluster_counts = out["scaffold_cluster"].value_counts()
    out["diversity_score"] = 1.0 / out["scaffold_cluster"].map(cluster_counts)
    out["pred_sweet_prob"] = np.clip(0.15 + 0.8 * out["sweet_score_soft"], 0.0, 1.0)
    out["acquisition_score"] = (
        0.5 * out["sweet_score_soft"]
        + 0.2 * out["pred_sweet_prob"]
        + 0.3 * out["diversity_score"]
    )

    return out.sort_values("acquisition_score", ascending=False).reset_index(drop=True)


def metric_card_row(cards: list[tuple[str, str]]) -> None:
    cols = st.columns(len(cards))
    for col, (label, value) in zip(cols, cards):
        col.metric(label, value)


def plot_pair_scatter(pair_df: pd.DataFrame) -> go.Figure:
    fig = px.scatter(
        pair_df,
        x="promiscuity",
        y="stereoselectivity_fraction",
        color="sweet_spot",
        hover_data=[c for c in ["pair_id", "scaffold_cluster"] if c in pair_df.columns],
        title="Promiscuity vs stereoselectivity fraction",
    )
    fig.add_vline(x=st.session_state.config["sweet_spot_min_sites"], line_dash="dash")
    fig.add_vline(x=st.session_state.config["sweet_spot_max_sites"], line_dash="dash")
    fig.add_hline(y=st.session_state.config["sweet_spot_min_fraction"], line_dash="dash")
    return fig


def plot_delta_auc_hist(hits_df: pd.DataFrame) -> go.Figure:
    return px.histogram(hits_df, x="delta_auc", nbins=40, title="ΔAUC distribution")


def plot_log2cr_hist(hits_df: pd.DataFrame) -> go.Figure:
    return px.histogram(hits_df, x="max_log2_cr", nbins=40, title="max log₂(CR) distribution")


def show_header() -> None:
    st.title(APP_NAME)
    st.caption(APP_TAGLINE)


def sidebar_controls() -> str:
    st.sidebar.header("Workflow")

    mode = st.sidebar.radio(
        "Mode",
        ["Analyze Real Data", "Simulation and Benchmarking", "Library Design"],
        index=0,
    )

    st.sidebar.divider()
    st.sidebar.subheader("Session")

    c1, c2 = st.sidebar.columns(2)
    if c1.button("Load demo", key="sidebar_load_demo"):
        st.session_state.raw_df = demo_real_dataset()
    if c2.button("Reset", key="sidebar_reset"):
        reset_session()
        st.rerun()

    st.sidebar.divider()
    st.sidebar.subheader("Global thresholds")

    cfg = st.session_state.config
    cfg["liganded_log2cr_threshold"] = st.sidebar.slider(
        "Liganded threshold, log₂(CR)",
        0.0, 4.0, float(cfg["liganded_log2cr_threshold"]), 0.1,
    )
    cfg["enantioselective_abs_delta_auc_threshold"] = st.sidebar.slider(
        "|ΔAUC| threshold",
        1.0, 50.0, float(cfg["enantioselective_abs_delta_auc_threshold"]), 1.0,
    )
    cfg["sweet_spot_min_sites"] = st.sidebar.number_input(
        "Sweet-spot min promiscuity",
        min_value=0,
        max_value=200,
        value=int(cfg["sweet_spot_min_sites"]),
    )
    cfg["sweet_spot_max_sites"] = st.sidebar.number_input(
        "Sweet-spot max promiscuity",
        min_value=1,
        max_value=1000,
        value=int(cfg["sweet_spot_max_sites"]),
    )
    cfg["sweet_spot_min_fraction"] = st.sidebar.slider(
        "Sweet-spot min stereoselectivity fraction",
        0.0, 1.0, float(cfg["sweet_spot_min_fraction"]), 0.01,
    )

    return mode

def page_home_real_data() -> None:
    metric_card_row(
        [
            ("Mode", "Analyze Real Data"),
            ("Template columns", str(len(REQUIRED_COLUMNS_REAL))),
            ("Core metrics", "CR, AUC, ΔAUC"),
        ]
    )

    st.markdown(
        "CysSelect is a browser-based research app for analyzing covalent fragment "
        "chemoproteomic screening data. It helps users identify liganded cysteines, "
        "compare enantiomer-specific behavior, prioritize higher-quality hits, explore "
        "simulated benchmark datasets, and rank virtual compounds for next-round library design."
    )

    st.info(
        "CysSelect is intended for research, method development, and educational use. "
        "Its outputs are computational summaries and prioritization aids, and should be "
        "interpreted alongside experimental validation and scientific judgment."
    )

    c1, c2, c3 = st.columns(3)
    c1.info(
        "**Upload and analyze screening data**\n\n"
        "Import chemoproteomic results, check data quality, and compute core metrics such as CR, AUC, and ΔAUC."
    )
    c2.info(
        "**Compare enantiomers and rank hits**\n\n"
        "Evaluate R and S enantiomer behavior, identify stereoselective interactions, and summarize pair-level promiscuity and selectivity."
    )
    c3.info(
        "**Simulate and design better libraries**\n\n"
        "Generate benchmark datasets, stress-test thresholds, and score virtual compound libraries to guide next-round selection."
    )

    st.subheader("Expected input")
    st.markdown(
        "CysSelect accepts tabular chemoproteomic screening data with required identifiers "
        "for compound pair, enantiomer, cysteine site, concentration, replicate, and control versus treatment intensities."
    )
    st.code(", ".join(REQUIRED_COLUMNS_REAL), language="text")

    st.subheader("Optional annotation fields")
    st.markdown(
        "Additional fields such as p-value, protein identifier, gene name, and scaffold cluster "
        "can improve interpretation, filtering, and prioritization."
    )
    st.code(", ".join(OPTIONAL_COLUMNS_REAL), language="text")

    st.subheader("What you can do in CysSelect")
    st.markdown(
        "- Upload and quality-check chemoproteomic screening tables\n"
        "- Call liganded and enantioselective hits\n"
        "- Summarize pair-level promiscuity and stereoselectivity\n"
        "- Identify sweet-spot fragment pairs\n"
        "- Explore simulated benchmark datasets\n"
        "- Score virtual libraries for next-round design"
    )

def page_upload_qc() -> None:
    st.header("Upload and QC")

    uploaded = st.file_uploader("Upload CSV", type=["csv"], key="upload_real_csv")
    if uploaded is not None:
        st.session_state.raw_df = pd.read_csv(uploaded)

    df = st.session_state.raw_df
    if df is None:
        st.info("Upload a CSV or click 'Load demo' in the sidebar.")
        return

    qc = validate_real_data(df)

    metric_card_row(
        [
            ("Rows", f"{qc['n_rows']:,}"),
            ("Unique pairs", f"{qc['n_pairs']:,}"),
            ("Unique cysteines", f"{qc['n_cysteines']:,}"),
            ("Missing fraction", f"{qc['missing_fraction']:.2%}"),
        ]
    )

    st.subheader("Preview")
    st.dataframe(df.head(20), use_container_width=True)

    if qc["errors"]:
        st.error("Dataset is not ready.")
        for e in qc["errors"]:
            st.write(f"• {e}")
    elif qc["warnings"]:
        st.warning("Dataset is usable with warnings.")
        for w in qc["warnings"]:
            st.write(f"• {w}")
    else:
        st.success("Dataset is ready for hit calling.")

    if qc["ready"]:
        safe_download_button(
            df,
            "Download current dataset",
            "cysselect_uploaded_data.csv",
            key="upload_qc_download_current_dataset",
        )


def page_hit_calling() -> None:
    st.header("Hit Calling")

    df = st.session_state.raw_df
    if df is None:
        st.info("No dataset loaded.")
        return

    if st.button("Run hit calling", type="primary", key="run_hit_calling"):
        hits_df, pair_summary_df = summarize_hits(
            df,
            liganded_log2cr_threshold=st.session_state.config["liganded_log2cr_threshold"],
            enantioselective_abs_delta_auc_threshold=st.session_state.config["enantioselective_abs_delta_auc_threshold"],
        )
        st.session_state.processed_hits_df = hits_df
        st.session_state.pair_summary_df = pair_summary_df

    hits_df = st.session_state.processed_hits_df
    if hits_df is None:
        st.info("Click 'Run hit calling' to compute CR, AUC, and ΔAUC.")
        return

    n_liganded = int(hits_df["liganded"].sum())
    n_enantio = int(hits_df["selectivity_class"].isin(["R_selective", "S_selective"]).sum())

    metric_card_row(
        [
            ("Interactions", f"{len(hits_df):,}"),
            ("Liganded", f"{n_liganded:,}"),
            ("Enantioselective", f"{n_enantio:,}"),
            ("% enantioselective among liganded", f"{(n_enantio / max(n_liganded, 1)):.1%}"),
        ]
    )

    c1, c2 = st.columns(2)
    with c1:
        st.plotly_chart(plot_log2cr_hist(hits_df), use_container_width=True)
    with c2:
        st.plotly_chart(plot_delta_auc_hist(hits_df), use_container_width=True)

    st.subheader("Processed hit table")
    st.dataframe(hits_df, use_container_width=True)
    safe_download_button(
        hits_df,
        "Download processed hits",
        "cysselect_hits.csv",
        key="hit_calling_download_processed_hits",
    )


def page_enantiomer_comparison() -> None:
    st.header("Enantiomer Comparison")

    hits_df = st.session_state.processed_hits_df
    pair_summary_df = st.session_state.pair_summary_df

    if hits_df is None or pair_summary_df is None:
        st.info("Run hit calling first.")
        return

    selected_pair = st.selectbox("Select pair", pair_summary_df["pair_id"].tolist(), key="enantiomer_pair_select")

    pair_meta = pair_summary_df[pair_summary_df["pair_id"] == selected_pair].iloc[0]
    pair_hits = hits_df[hits_df["pair_id"] == selected_pair].copy()

    metric_card_row(
        [
            ("Promiscuity", str(int(pair_meta["promiscuity"]))),
            ("Stereo fraction", f"{pair_meta['stereoselectivity_fraction']:.2f}"),
            ("R-selective", str(int(pair_meta["n_r_selective"]))),
            ("S-selective", str(int(pair_meta["n_s_selective"]))),
        ]
    )

    st.write(f"Sweet spot: **{bool(pair_meta['sweet_spot'])}**")

    fig = px.scatter(
        pair_hits,
        x="auc_r",
        y="auc_s",
        color="selectivity_class",
        hover_data=[c for c in ["cysteine_id", "gene_name", "max_log2_cr", "delta_auc"] if c in pair_hits.columns],
        title=f"R vs S AUC for {selected_pair}",
    )
    max_x = max(float(pair_hits["auc_r"].fillna(0).max()), 1.0)
    max_y = max(float(pair_hits["auc_s"].fillna(0).max()), 1.0)
    lim = max(max_x, max_y)
    fig.add_shape(type="line", x0=0, y0=0, x1=lim, y1=lim)
    st.plotly_chart(fig, use_container_width=True)

    top_hits = pair_hits.sort_values("delta_auc", key=np.abs, ascending=False)
    st.subheader("Top differential cysteines")
    st.dataframe(top_hits, use_container_width=True)

    safe_download_button(
        top_hits,
        "Download pair details",
        f"{selected_pair}_details.csv",
        key=f"enantiomer_download_{selected_pair}",
    )


def page_pair_summary() -> None:
    st.header("Pair-Level Summary")

    pair_summary_df = st.session_state.pair_summary_df
    if pair_summary_df is None:
        st.info("Run hit calling first.")
        return

    metric_card_row(
        [
            ("Pairs", f"{len(pair_summary_df):,}"),
            ("Median promiscuity", f"{pair_summary_df['promiscuity'].median():.1f}"),
            ("Median stereo fraction", f"{pair_summary_df['stereoselectivity_fraction'].median():.2f}"),
            ("Sweet-spot pairs", f"{int(pair_summary_df['sweet_spot'].sum())}"),
        ]
    )

    st.plotly_chart(plot_pair_scatter(pair_summary_df), use_container_width=True)

    c1, c2 = st.columns(2)
    with c1:
        st.plotly_chart(
            px.histogram(pair_summary_df, x="promiscuity", nbins=30, title="Promiscuity distribution"),
            use_container_width=True,
        )
    with c2:
        st.plotly_chart(
            px.histogram(pair_summary_df, x="stereoselectivity_fraction", nbins=30, title="Stereoselectivity fraction distribution"),
            use_container_width=True,
        )

    st.subheader("Ranked pairs")
    ranked = pair_summary_df.sort_values(
        ["sweet_spot", "stereoselectivity_fraction", "promiscuity"],
        ascending=[False, False, True],
    )
    st.dataframe(ranked, use_container_width=True)

    safe_download_button(
        ranked,
        "Download pair summary",
        "cysselect_pair_summary.csv",
        key="pair_summary_download",
    )


def page_hit_prioritization() -> None:
    st.header("Hit Prioritization")

    pair_summary_df = st.session_state.pair_summary_df
    if pair_summary_df is None:
        st.info("Run hit calling first.")
        return

    st.subheader("Weights")
    c1, c2, c3, c4 = st.columns(4)
    w_prom = c1.slider("Promiscuity desirability", 0.0, 1.0, 0.25, 0.05, key="w_prom")
    w_stereo = c2.slider("Stereoselectivity reward", 0.0, 1.0, 0.35, 0.05, key="w_stereo")
    w_sweet = c3.slider("Sweet-spot reward", 0.0, 1.0, 0.25, 0.05, key="w_sweet")
    w_div = c4.slider("Diversity reward", 0.0, 1.0, 0.15, 0.05, key="w_div")

    ranked = pair_summary_df.copy()

    def smooth_range_score(x, low=4, high=64, width=8):
        left = 1 / (1 + np.exp(-(x - low) / width))
        right = 1 / (1 + np.exp((x - high) / width))
        return left * right

    ranked["prom_score"] = smooth_range_score(
        ranked["promiscuity"],
        st.session_state.config["sweet_spot_min_sites"],
        st.session_state.config["sweet_spot_max_sites"],
        8,
    )

    if "scaffold_cluster" in ranked.columns:
        cluster_counts = ranked["scaffold_cluster"].fillna(-1).value_counts()
        ranked["diversity_score"] = 1.0 / ranked["scaffold_cluster"].fillna(-1).map(cluster_counts)
    else:
        ranked["diversity_score"] = 1.0

    ranked["priority_score"] = (
        w_prom * ranked["prom_score"]
        + w_stereo * ranked["stereoselectivity_fraction"]
        + w_sweet * ranked["sweet_spot"].astype(float)
        + w_div * ranked["diversity_score"]
    )

    ranked = ranked.sort_values("priority_score", ascending=False).reset_index(drop=True)

    metric_card_row(
        [
            ("Top pair", str(ranked.iloc[0]["pair_id"])),
            ("Top score", f"{ranked.iloc[0]['priority_score']:.3f}"),
            ("High-confidence sweet-spot pairs", f"{int(ranked['sweet_spot'].sum())}"),
            ("Mean priority score", f"{ranked['priority_score'].mean():.3f}"),
        ]
    )

    st.plotly_chart(
        px.bar(ranked.head(20), x="pair_id", y="priority_score", title="Top prioritized pairs"),
        use_container_width=True,
    )

    st.dataframe(ranked, use_container_width=True)

    safe_download_button(
        ranked,
        "Download prioritized hits",
        "cysselect_prioritized_pairs.csv",
        key="prioritization_download",
    )


def page_simulation_setup() -> None:
    st.header("Simulation and Benchmarking")

    c1, c2, c3 = st.columns(3)
    n_pairs = c1.number_input("Number of enantiopairs", min_value=10, max_value=500, value=DEFAULT_CONFIG["n_sim_pairs"], key="sim_n_pairs")
    n_cys = c2.number_input("Number of cysteines", min_value=50, max_value=10000, value=DEFAULT_CONFIG["n_sim_cysteines"], key="sim_n_cys")
    seed = c3.number_input("Seed", min_value=1, max_value=999999, value=DEFAULT_CONFIG["seed"], key="sim_seed")

    c4, c5 = st.columns(2)
    noise_sd = c4.slider("Assay noise SD", 0.01, 1.0, DEFAULT_CONFIG["sim_noise_sd"], 0.01, key="sim_noise")
    missing_rate = c5.slider("Missing-data rate", 0.0, 0.5, DEFAULT_CONFIG["sim_missing_rate"], 0.01, key="sim_missing")

    if st.button("Generate simulated dataset", type="primary", key="generate_sim_dataset"):
        sim_df = simulate_dataset(
            n_pairs=int(n_pairs),
            n_cysteines=int(n_cys),
            noise_sd=float(noise_sd),
            missing_rate=float(missing_rate),
            seed=int(seed),
        )
        st.session_state.sim_raw_df = sim_df

        sim_hits_df, sim_pair_summary_df = summarize_hits(
            sim_df,
            liganded_log2cr_threshold=st.session_state.config["liganded_log2cr_threshold"],
            enantioselective_abs_delta_auc_threshold=st.session_state.config["enantioselective_abs_delta_auc_threshold"],
        )
        st.session_state.sim_hits_df = sim_hits_df
        st.session_state.sim_pair_summary_df = sim_pair_summary_df

    sim_df = st.session_state.sim_raw_df
    if sim_df is None:
        st.info("Generate a simulated dataset to continue.")
        return

    metric_card_row(
        [
            ("Sim rows", f"{len(sim_df):,}"),
            ("Sim pairs", f"{sim_df['pair_id'].nunique():,}"),
            ("Sim cysteines", f"{sim_df['cysteine_id'].nunique():,}"),
            ("Concentrations", ", ".join(map(str, sorted(sim_df["concentration_uM"].unique())))),
        ]
    )

    st.subheader("Simulated raw preview")
    st.dataframe(sim_df.head(20), use_container_width=True)

    sim_hits_df = st.session_state.sim_hits_df
    sim_pair_summary_df = st.session_state.sim_pair_summary_df

    if sim_hits_df is not None and sim_pair_summary_df is not None:
        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(plot_log2cr_hist(sim_hits_df), use_container_width=True)
        with c2:
            st.plotly_chart(plot_delta_auc_hist(sim_hits_df), use_container_width=True)

        st.plotly_chart(plot_pair_scatter(sim_pair_summary_df), use_container_width=True)

        safe_download_button(sim_df, "Download simulated raw data", "cysselect_sim_raw.csv", key="sim_download_raw")
        safe_download_button(sim_hits_df, "Download simulated hits", "cysselect_sim_hits.csv", key="sim_download_hits")
        safe_download_button(sim_pair_summary_df, "Download simulated pair summary", "cysselect_sim_pair_summary.csv", key="sim_download_pair_summary")


def page_library_design() -> None:
    st.header("Library Design")

    c1, c2 = st.columns(2)
    n_candidates = c1.number_input("Number of virtual candidates", min_value=50, max_value=5000, value=500, key="lib_n_candidates")
    seed = c2.number_input("Seed for virtual library", min_value=1, max_value=999999, value=42, key="lib_seed")

    if st.button("Generate virtual library", key="generate_virtual_library"):
        st.session_state.virtual_library_df = generate_virtual_library(n_candidates=int(n_candidates), seed=int(seed))

    uploaded = st.file_uploader("Or upload virtual library CSV", type=["csv"], key="virtual_upload")
    if uploaded is not None:
        st.session_state.virtual_library_df = pd.read_csv(uploaded)

    virtual_df = st.session_state.virtual_library_df
    if virtual_df is None:
        st.info("Generate or upload a virtual library.")
        return

    st.subheader("Virtual library preview")
    st.dataframe(virtual_df.head(20), use_container_width=True)

    if st.button("Score virtual library", type="primary", key="score_virtual_library"):
        st.session_state.scored_library_df = score_virtual_library(virtual_df)

    scored_df = st.session_state.scored_library_df
    if scored_df is None:
        return

    metric_card_row(
        [
            ("Candidates scored", f"{len(scored_df):,}"),
            ("Best acquisition score", f"{scored_df['acquisition_score'].max():.3f}"),
            ("Median predicted promiscuity", f"{scored_df['pred_promiscuity'].median():.1f}"),
            ("Median predicted stereo fraction", f"{scored_df['pred_stereo_fraction'].median():.2f}"),
        ]
    )

    fig = px.scatter(
        scored_df,
        x="pred_promiscuity",
        y="pred_stereo_fraction",
        color="acquisition_score",
        hover_data=["candidate_id", "scaffold_cluster", "sweet_score_soft"],
        title="Virtual library scoring",
    )
    fig.add_vline(x=st.session_state.config["sweet_spot_min_sites"], line_dash="dash")
    fig.add_vline(x=st.session_state.config["sweet_spot_max_sites"], line_dash="dash")
    fig.add_hline(y=st.session_state.config["sweet_spot_min_fraction"], line_dash="dash")
    st.plotly_chart(fig, use_container_width=True)

    selected_n = st.slider("Number to select", 5, min(100, len(scored_df)), 38, 1, key="selected_n")
    selected = scored_df.head(selected_n).copy()

    st.subheader("Selected next-round batch")
    st.dataframe(selected, use_container_width=True)

    safe_download_button(
        scored_df,
        "Download scored library",
        "cysselect_scored_library.csv",
        key="library_design_download_scored",
    )
    safe_download_button(
        selected,
        "Download selected batch",
        "cysselect_selected_batch.csv",
        key="library_design_download_selected",
    )


def page_export() -> None:
    st.header("Reports and Export")

    tabs = st.tabs(["Raw", "Processed hits", "Pair summary", "Simulation", "Library design"])

    with tabs[0]:
        if st.session_state.raw_df is not None:
            safe_download_button(
                st.session_state.raw_df,
                "Download raw uploaded data",
                "cysselect_raw.csv",
                key="export_raw",
            )
            st.dataframe(st.session_state.raw_df.head(20), use_container_width=True)
        else:
            st.info("No raw uploaded data.")

    with tabs[1]:
        if st.session_state.processed_hits_df is not None:
            safe_download_button(
                st.session_state.processed_hits_df,
                "Download processed hits",
                "cysselect_hits.csv",
                key="export_hits",
            )
            st.dataframe(st.session_state.processed_hits_df.head(20), use_container_width=True)
        else:
            st.info("No processed hits yet.")

    with tabs[2]:
        if st.session_state.pair_summary_df is not None:
            safe_download_button(
                st.session_state.pair_summary_df,
                "Download pair summary",
                "cysselect_pair_summary.csv",
                key="export_pair_summary",
            )
            st.dataframe(st.session_state.pair_summary_df.head(20), use_container_width=True)
        else:
            st.info("No pair summary yet.")

    with tabs[3]:
        if st.session_state.sim_pair_summary_df is not None:
            safe_download_button(
                st.session_state.sim_pair_summary_df,
                "Download simulation summary",
                "cysselect_sim_summary.csv",
                key="export_sim_summary",
            )
            st.dataframe(st.session_state.sim_pair_summary_df.head(20), use_container_width=True)
        else:
            st.info("No simulation summary yet.")

    with tabs[4]:
        if st.session_state.scored_library_df is not None:
            safe_download_button(
                st.session_state.scored_library_df,
                "Download scored library",
                "cysselect_scored_library.csv",
                key="export_scored_library",
            )
            st.dataframe(st.session_state.scored_library_df.head(20), use_container_width=True)
        else:
            st.info("No scored library yet.")

    st.subheader("Session configuration")
    st.json(st.session_state.config)


def run_app() -> None:
    init_state()
    show_header()
    mode = sidebar_controls()

    if mode == "Analyze Real Data":
        tabs = st.tabs(
            [
                "Home",
                "Upload and QC",
                "Hit Calling",
                "Enantiomer Comparison",
                "Pair-Level Summary",
                "Hit Prioritization",
                "Export",
            ]
        )
        with tabs[0]:
            page_home_real_data()
        with tabs[1]:
            page_upload_qc()
        with tabs[2]:
            page_hit_calling()
        with tabs[3]:
            page_enantiomer_comparison()
        with tabs[4]:
            page_pair_summary()
        with tabs[5]:
            page_hit_prioritization()
        with tabs[6]:
            page_export()

    elif mode == "Simulation and Benchmarking":
        tabs = st.tabs(["Simulation Setup", "Export"])
        with tabs[0]:
            page_simulation_setup()
        with tabs[1]:
            page_export()

    elif mode == "Library Design":
        tabs = st.tabs(["Library Design", "Export"])
        with tabs[0]:
            page_library_design()
        with tabs[1]:
            page_export()


if __name__ == "__main__":
    run_app()
