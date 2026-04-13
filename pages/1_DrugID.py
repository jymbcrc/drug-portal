from pathlib import Path
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# =========================
# Page setup
# =========================
st.set_page_config(
    page_title="Drug Proteomics Portal",
    layout="wide"
)

# =========================
# Paths
# =========================
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "drugs"

# =========================
# Helpers
# =========================
@st.cache_data
def get_drug_list(data_dir: Path):
    if not data_dir.exists():
        return []
    return sorted([p.name for p in data_dir.iterdir() if p.is_dir()])


@st.cache_data
def load_summary_json(drug_dir_str: str):
    drug_dir = Path(drug_dir_str)
    file = drug_dir / "summary.json"
    if not file.exists():
        return {}
    with open(file, "r", encoding="utf-8") as f:
        return json.load(f)


@st.cache_data
def load_csv(file_path_str: str):
    file_path = Path(file_path_str)
    if not file_path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(file_path)
    except Exception:
        return pd.DataFrame()


@st.cache_data
def load_markdown(file_path_str: str):
    file_path = Path(file_path_str)
    if not file_path.exists():
        return ""
    return file_path.read_text(encoding="utf-8")


def show_metric_with_caption(label, value, caption):
    st.metric(label, value)
    if caption:
        st.caption(caption)


def plot_volcano_interactive(
    df,
    drug_name,
    fc_thresh=0.5,
    fdr_thresh=0.05,
    top_n_labels=12,
    highlight_gene=None
):
    if df.empty:
        return None, pd.DataFrame()

    sub = df.copy()

    if "drug" in sub.columns:
        sub = sub.loc[sub["drug"] == drug_name].copy()

    required_cols = ["protein", "log2FC", "p_value", "p_value_adj"]
    for c in required_cols:
        if c not in sub.columns:
            return None, pd.DataFrame()

    sub = sub.dropna(subset=["log2FC", "p_value_adj"])

    if sub.empty:
        return None, pd.DataFrame()

    if "gene_symbol" in sub.columns:
        sub["label_name"] = sub["gene_symbol"].fillna("").astype(str)
        sub.loc[sub["label_name"].str.strip() == "", "label_name"] = sub["protein"].astype(str)
    else:
        sub["label_name"] = sub["protein"].astype(str)

    sub["p_value_adj"] = pd.to_numeric(sub["p_value_adj"], errors="coerce")
    sub["p_value"] = pd.to_numeric(sub["p_value"], errors="coerce")
    sub["log2FC"] = pd.to_numeric(sub["log2FC"], errors="coerce")

    if "t_stat" in sub.columns:
        sub["t_stat"] = pd.to_numeric(sub["t_stat"], errors="coerce")
    if "mean_drug" in sub.columns:
        sub["mean_drug"] = pd.to_numeric(sub["mean_drug"], errors="coerce")
    if "mean_DMSO" in sub.columns:
        sub["mean_DMSO"] = pd.to_numeric(sub["mean_DMSO"], errors="coerce")

    sub = sub.dropna(subset=["log2FC", "p_value_adj", "p_value"])

    if sub.empty:
        return None, pd.DataFrame()

    sub["p_value_adj"] = sub["p_value_adj"].clip(lower=1e-300)
    sub["p_value"] = sub["p_value"].clip(lower=1e-300)
    sub["neglog10P"] = -np.log10(sub["p_value_adj"])

    sub["sig_group"] = "Not significant"
    sub.loc[
        (sub["log2FC"] >= fc_thresh) & (sub["p_value_adj"] < fdr_thresh),
        "sig_group"
    ] = "Up"
    sub.loc[
        (sub["log2FC"] <= -fc_thresh) & (sub["p_value_adj"] < fdr_thresh),
        "sig_group"
    ] = "Down"

    sig_for_label = sub.loc[sub["sig_group"] != "Not significant"].copy()
    sig_for_label["label_score"] = sig_for_label["neglog10P"] * sig_for_label["log2FC"].abs()
    label_df = sig_for_label.nlargest(top_n_labels, "label_score")

    color_map = {
        "Not significant": "#BDBDBD",
        "Up": "#E64B35",
        "Down": "#3C91E6"
    }

    fig = go.Figure()

    for group in ["Not significant", "Up", "Down"]:
        plot_df = sub[sub["sig_group"] == group].copy()

        hover_text = (
            "Gene: " + plot_df["label_name"].astype(str) +
            "<br>Protein: " + plot_df["protein"].astype(str) +
            "<br>log2FC: " + plot_df["log2FC"].round(3).astype(str) +
            "<br>adj.P: " + plot_df["p_value_adj"].apply(lambda x: f"{x:.2e}") +
            "<br>P-value: " + plot_df["p_value"].apply(lambda x: f"{x:.2e}")
        )

        if "t_stat" in plot_df.columns:
            hover_text += "<br>t-stat: " + plot_df["t_stat"].round(3).astype(str)

        if "mean_drug" in plot_df.columns:
            hover_text += "<br>mean_drug: " + plot_df["mean_drug"].round(3).astype(str)

        if "mean_DMSO" in plot_df.columns:
            hover_text += "<br>mean_DMSO: " + plot_df["mean_DMSO"].round(3).astype(str)

        fig.add_trace(go.Scattergl(
            x=plot_df["log2FC"],
            y=plot_df["neglog10P"],
            mode="markers",
            name=group,
            text=plot_df["label_name"],
            hovertext=hover_text,
            hoverinfo="text",
            marker=dict(
                color=color_map[group],
                size=7,
                opacity=0.75
            )
        ))

    if len(label_df) > 0:
        fig.add_trace(go.Scatter(
            x=label_df["log2FC"],
            y=label_df["neglog10P"],
            mode="text",
            text=label_df["label_name"],
            textposition="top center",
            showlegend=False,
            name="Labels"
        ))

    if highlight_gene is not None and str(highlight_gene).strip() != "":
        highlight_df = sub[
            sub["label_name"].str.upper() == highlight_gene.strip().upper()
        ].copy()

        if not highlight_df.empty:
            fig.add_trace(go.Scatter(
                x=highlight_df["log2FC"],
                y=highlight_df["neglog10P"],
                mode="markers+text",
                text=highlight_df["label_name"],
                textposition="top center",
                name="Highlighted gene",
                marker=dict(
                    size=14,
                    symbol="circle-open",
                    line=dict(width=2, color="black"),
                    color="gold"
                ),
                hoverinfo="skip"
            ))

    fig.add_vline(x=fc_thresh, line_dash="dash", line_width=1)
    fig.add_vline(x=-fc_thresh, line_dash="dash", line_width=1)
    fig.add_hline(y=-np.log10(fdr_thresh), line_dash="dash", line_width=1)

    fig.update_layout(
        title=f"Volcano plot: {drug_name}",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(adjusted P-value)",
        template="simple_white",
        height=800,
        legend_title="Group",
        xaxis=dict(
            title_font=dict(size=18),
            tickfont=dict(size=18)
        ),  
        yaxis=dict(
            title_font=dict(size=18),
            tickfont=dict(size=18)
        ),
        legend=dict(
            title_font=dict(size=18),
            font=dict(size=18)
        )
    )

    return fig, sub


# =========================
# Main title
# =========================
st.title("Drug Proteomics Portal")
st.write("Interactive companion portal for browsing drug-level proteomic perturbation results.")

# =========================
# Check data dir
# =========================
if not DATA_DIR.exists():
    st.error(f"Cannot find data directory: {DATA_DIR}")
    st.stop()

drug_list = get_drug_list(DATA_DIR)

if len(drug_list) == 0:
    st.warning("No drug folders were found under data/drugs/")
    st.stop()

# =========================
# Sidebar
# =========================
st.sidebar.header("Browse")

search_text = st.sidebar.text_input("Search drug name", "")

filtered_drugs = [
    d for d in drug_list
    if search_text.lower() in d.lower()
]

if len(filtered_drugs) == 0:
    st.sidebar.warning("No matching drugs found.")
    st.stop()

selected_drug = st.sidebar.selectbox("Select a drug", filtered_drugs)

drug_dir = DATA_DIR / selected_drug

summary = load_summary_json(str(drug_dir))
diff_df = load_csv(str(drug_dir / "diff_proteins.csv"))
hallmark_df = load_csv(str(drug_dir / "hallmark.csv"))
go_df = load_csv(str(drug_dir / "go_bp.csv"))
ttest_df = load_csv(str(drug_dir / "ttest.csv"))
ai_summary_md = load_markdown(str(drug_dir / "ai_summary.md"))

# =========================
# Header
# =========================
drug_name = summary.get("drug_name", selected_drug)
target = summary.get("target", "NA")
pathway = summary.get("pathway", "NA")
research_area = summary.get("research_area", "NA")

st.header(drug_name)

info_col1, info_col2, info_col3 = st.columns(3)
with info_col1:
    st.markdown(f"**Target:** {target}")
with info_col2:
    st.markdown(f"**Pathway:** {pathway}")
with info_col3:
    st.markdown(f"**Research area:** {research_area}")

# =========================
# Summary metrics
# =========================
st.subheader("Summary")

m1, m2, m3 = st.columns(3)

protein_criteria = summary.get("n_sig_proteins_criteria", {})
go_criteria = summary.get("n_sig_go_bp_criteria", {})
hallmark_criteria = summary.get("n_sig_hallmark_criteria", {})

protein_caption = []
if "FDR" in protein_criteria:
    protein_caption.append(f"FDR {protein_criteria['FDR']}")
if "log2FC" in protein_criteria and protein_criteria["log2FC"] != "None":
    protein_caption.append(f"|log2FC| {protein_criteria['log2FC']}")
protein_caption = ", ".join(protein_caption)

go_caption = []
if "FDR" in go_criteria:
    go_caption.append(f"FDR {go_criteria['FDR']}")
go_caption = ", ".join(go_caption)

hallmark_caption = []
if "FDR" in hallmark_criteria and hallmark_criteria["FDR"] != "None":
    hallmark_caption.append(f"FDR {hallmark_criteria['FDR']}")
elif "FDR" in hallmark_criteria and hallmark_criteria["FDR"] == "None":
    hallmark_caption.append("No FDR filter")
hallmark_caption = ", ".join(hallmark_caption)

with m1:
    show_metric_with_caption(
        "Significant proteins",
        summary.get("n_sig_proteins_total", "NA"),
        protein_caption
    )

with m2:
    show_metric_with_caption(
        "Significant GO BP pathways",
        summary.get("n_sig_go_bp_total", "NA"),
        go_caption
    )

with m3:
    show_metric_with_caption(
        "Significant Hallmark pathways",
        summary.get("n_sig_hallmark_total", "NA"),
        hallmark_caption
    )

# =========================
# Tabs
# =========================
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Volcano plot",
    "ID card",
    "All protein statistics",
    "Hallmark",
    "GO BP"
])

# =========================
# Tab 1: Volcano plot
# =========================
with tab1:
    if ttest_df.empty:
        st.info("No ttest.csv found.")
    else:
        v1, v2, v3, v4 = st.columns([1, 1, 1, 1.2])

        with v1:
            fc_thresh = st.slider(
                "log2FC threshold",
                min_value=0.0,
                max_value=2.0,
                value=0.5,
                step=0.1,
                key="volcano_fc"
            )

        with v2:
            fdr_thresh = st.select_slider(
                "Adjusted P threshold",
                options=[0.1, 0.05, 0.01, 0.001],
                value=0.05,
                key="volcano_fdr"
            )

        with v3:
            top_n_labels = st.slider(
                "Top labels",
                min_value=0,
                max_value=30,
                value=12,
                step=1,
                key="volcano_labels"
            )

        with v4:
            highlight_gene = st.text_input(
                "Highlight gene",
                value="",
                key="volcano_highlight"
            )

        fig, volcano_df = plot_volcano_interactive(
            df=ttest_df,
            drug_name=selected_drug,
            fc_thresh=fc_thresh,
            fdr_thresh=fdr_thresh,
            top_n_labels=top_n_labels,
            highlight_gene=highlight_gene
        )

        if fig is None or volcano_df.empty:
            st.warning("No valid volcano plot data available.")
        else:
            st.plotly_chart(fig, use_container_width=True)

            n_up = ((volcano_df["log2FC"] >= fc_thresh) & (volcano_df["p_value_adj"] < fdr_thresh)).sum()
            n_down = ((volcano_df["log2FC"] <= -fc_thresh) & (volcano_df["p_value_adj"] < fdr_thresh)).sum()
            n_total = volcano_df.shape[0]

            s1, s2, s3 = st.columns(3)
            s1.metric("Total proteins", int(n_total))
            s2.metric("Upregulated", int(n_up))
            s3.metric("Downregulated", int(n_down))

            st.markdown("### Significant proteins / genes")

            sig_cols = ["sig_group", "label_name", "protein", "log2FC", "p_value_adj", "p_value"]
            extra_cols = [c for c in ["t_stat", "mean_drug", "mean_DMSO"] if c in volcano_df.columns]
            sig_cols = sig_cols + extra_cols

            sig_df = volcano_df.loc[
                volcano_df["sig_group"] != "Not significant",
                sig_cols
            ].sort_values(["p_value_adj", "log2FC"], ascending=[True, False])

            st.write(f"Rows shown: {len(sig_df)}")
            st.dataframe(sig_df, use_container_width=True)

            csv_data = sig_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download displayed volcano significant table CSV",
                data=csv_data,
                file_name=f"{selected_drug}_volcano_sig.csv",
                mime="text/csv"
            )

# =========================
# Tab 2: Summary
# =========================
with tab2:
    if ai_summary_md.strip():
        st.markdown(ai_summary_md)
    else:
        st.info("No ai_summary.md found.")

# =========================
# Tab 3: Identified proteins
# =========================
with tab3:
    st.subheader("Identified proteins")

    if diff_df.empty:
        st.info("No diff_proteins.csv found.")
    else:
        if "direction" in diff_df.columns:
            dir_option = st.radio(
                "Filter by direction",
                ["All", "up", "down"],
                horizontal=True
            )

            show_df = diff_df.copy()
            if dir_option != "All":
                show_df = show_df.loc[show_df["direction"] == dir_option].copy()
        else:
            show_df = diff_df.copy()

        st.write(f"Rows shown: {len(show_df)}")
        st.dataframe(show_df, use_container_width=True)

        csv_data = show_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download displayed differential proteins CSV",
            data=csv_data,
            file_name=f"{selected_drug}_diff_proteins.csv",
            mime="text/csv"
        )

# =========================
# Tab 4: Hallmark
# =========================
with tab4:
    st.subheader("Hallmark pathway profile")

    if hallmark_df.empty:
        st.info("No hallmark.csv found.")
    else:
        st.dataframe(hallmark_df, use_container_width=True)

        csv_data = hallmark_df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download Hallmark CSV",
            data=csv_data,
            file_name=f"{selected_drug}_hallmark.csv",
            mime="text/csv"
        )

# =========================
# Tab 5: GO BP
# =========================
with tab5:
    st.subheader("GO Biological Process pathways")

    if go_df.empty:
        st.info("No go_bp.csv found.")
    else:
        if "direction" in go_df.columns:
            go_dir = st.radio(
                "Filter GO direction",
                ["All", "up", "down"],
                horizontal=True,
                key="go_direction"
            )

            show_go = go_df.copy()
            if go_dir != "All":
                show_go = show_go.loc[show_go["direction"] == go_dir].copy()
        else:
            show_go = go_df.copy()

        st.write(f"Rows shown: {len(show_go)}")
        st.dataframe(show_go, use_container_width=True)

        csv_data = show_go.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download displayed GO BP CSV",
            data=csv_data,
            file_name=f"{selected_drug}_go_bp.csv",
            mime="text/csv"
        )

# =========================
# Footer
# =========================
st.divider()
st.caption(f"Current folder: {drug_dir}")