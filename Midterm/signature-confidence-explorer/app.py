"""Signature Confidence Explorer: how far can you trust an OCTAD 2.0 disease signature?

Run locally:  streamlit run app.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

DATA = Path(__file__).parent / "data" / "processed"
CAP = 50  # -log10(adjusted p) is capped here for display; raw values are untouched
TIER_ORDER = ["low", "medium", "high"]
TIER_LABEL = {"low": "Low (1 dataset)", "medium": "Medium (2 to 5 datasets)", "high": "High (6 or more datasets)"}
# Three points taken from the Viridis scale so tiers match the colour of the k_datasets scale used in the volcano plot.
TIER_COLOR = {"low": "#440154", "medium": "#21918c", "high": "#7ad151"}
DEFAULT_DISEASE = "Diabetic-Nephropathies"

st.set_page_config(page_title="Signature Confidence Explorer", page_icon=":dna:", layout="wide")


# ----------------------------------------------------------------------------- data
@st.cache_data(show_spinner="Loading data...")
def load():
    ds = pd.read_csv(DATA / "disease_summary.csv")
    tested = ds["n_genes"] * (1 - ds["pct_genes_untested"] / 100)
    ds["pct_sig_of_tested"] = ds["n_sig"] / tested * 100
    ds["tier_label"] = ds["tier"].map(TIER_LABEL)
    tg = pd.read_csv(DATA / "top_genes.csv")
    tg["neglog10_padj"] = (-np.log10(tg["padj_meta"].clip(lower=1e-300))).clip(upper=CAP)
    tg["at_cap"] = tg["neglog10_padj"] >= CAP
    audit = json.loads((DATA / "audit.json").read_text())
    gs = pd.read_csv(DATA / "gene_sample.csv", usecols=["log2FoldChange_meta", "padj_meta", "k_datasets"])
    gs["abs_lfc"] = gs["log2FoldChange_meta"].abs()
    gs["tier"] = pd.cut(gs["k_datasets"], [0, 1, 5, np.inf], labels=TIER_ORDER).astype(str)
    return ds, tg, audit, gs


def pretty(name: str) -> str:
    return name.replace("-", " ")


def i2_label(v) -> str:
    """Higgins et al. (2003) rule of thumb: roughly 25 / 50 / 75 percent."""
    if pd.isna(v):
        return "Not measurable"
    return "Low" if v < 25 else "Moderate" if v < 50 else "Substantial" if v < 75 else "Considerable"


def fmt_p(p: float) -> str:
    return "0 (below display limit)" if p == 0 else f"{p:.2e}"


diseases, top_genes, audit, gene_sample = load()
n_dis = len(diseases)
tier_counts = diseases["tier"].value_counts().reindex(TIER_ORDER)
pct_low = tier_counts["low"] / n_dis * 100

# ----------------------------------------------------------------------------- header
st.title("Signature Confidence Explorer")
st.markdown("**How many datasets stand behind a disease signature, and do they agree?**")

tab_home, tab_disease, tab_compare, tab_data, tab_limits = st.tabs(
    ["Start here", "Explore a disease", "Compare diseases", "Data and preparation", "Limitations"]
)

# ----------------------------------------------------------------------------- tab 1
with tab_home:
    st.markdown("**For researchers about to build on a disease signature.** Is it backed by more than one study, and do the studies agree?")
    k1, k2, k3 = st.columns(3)
    k1.metric("Diseases on a single dataset", f"{pct_low:.0f}%")
    k1.caption(f"{tier_counts['low']} of {n_dis}")
    k2.metric("Diseases with 6+ datasets", f"{tier_counts['high'] / n_dis * 100:.0f}%")
    k2.caption(f"{tier_counts['high']} of {n_dis}")
    k3.metric("Gene results analysed", f"{audit['rows'] / 1e6:.1f} million")
    k3.caption(f"{n_dis} diseases")
    st.info("A disease menu makes every disease look equally solid. This app shows which ones are.")
    st.markdown(
        "**Three questions it answers**\n"
        "1. How much evidence stands behind each disease? *(Explore a disease, Compare diseases)*\n"
        "2. Do more datasets mean stronger or larger signals? *(Compare diseases)*\n"
        "3. Do more datasets bring more agreement? *(Compare diseases)*"
    )
    st.markdown("**Start:** open *Explore a disease*, pick a disease, read its confidence card.")
    with st.expander("What do the terms mean?"):
        st.markdown(
            "- **log2FC:** size and direction of change (+1 is about double).\n"
            "- **padj:** significance after correcting for many genes. Smaller is stronger.\n"
            "- **Datasets (k):** how many datasets back a gene's result.\n"
            "- **I²:** how much those datasets disagree (0 to 100%). Needs 2+ datasets."
        )
    st.caption("Data: OCTAD 2.0 random effects meta analysis. Descriptive only; no drug recommendations.")

# ----------------------------------------------------------------------------- tab 2
with tab_disease:
    options = sorted(diseases["disease"])
    c1, c2, c3 = st.columns([2, 1.5, 1.5])
    with c1:
        disease = st.selectbox(
            "Disease (click, then type to search)", options, index=options.index(DEFAULT_DISEASE), format_func=pretty,
            help="All 350 OCTAD 2.0 diseases.",
        )
    row = diseases.loc[diseases["disease"] == disease].iloc[0]
    sub = top_genes[top_genes["disease"] == disease]
    k_best = int(row["k_max"])
    multi = k_best >= 2
    with c2:
        min_k = st.slider(
            "Minimum datasets per gene", 1, max(k_best, 2) if multi else 2, 1, disabled=not multi,
            help="Hide genes supported by fewer datasets. Disabled when every gene has a single dataset.",
        )
    with c3:
        show = st.radio("Direction", ["Both", "Up in disease", "Down in disease"], horizontal=True)

    # Confidence card
    st.markdown(f"#### {pretty(disease)}")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Confidence tier", row["tier"].capitalize(), help="Based on the best supported gene: low = 1 dataset, medium = 2 to 5, high = 6 or more.")
    m1.caption(TIER_LABEL[row["tier"]])
    m2.metric("Best supported gene", f"{k_best} dataset{'s' if k_best != 1 else ''}")
    m2.caption(f"A typical gene has {row['k_median']:.0f}")
    m3.metric("Significant genes (padj < 0.05)", f"{int(row['n_sig']):,}")
    m3.caption(f"{row['pct_sig_of_tested']:.0f}% of the genes that were tested")
    m4.metric("Dataset disagreement (I²)", i2_label(row["median_I2_k2plus"]),
              help="Median I² of this disease's multi dataset genes. Low under 25%, moderate 25 to 50%, substantial 50 to 75%, considerable above 75%.")
    m4.caption(f"Median {row['median_I2_k2plus']:.0f}%" if multi else "Needs 2 or more datasets")

    if row["tier"] == "low":
        st.warning("One dataset behind every gene. Nothing has replicated it, and agreement can't be measured.")
    elif row["tier"] == "medium":
        st.info(f"Some replication (best gene: {k_best} datasets), but only {row['pct_genes_k6plus']:.0f}% of genes have 6+. Check your genes.")
    else:
        st.success(f"Well supported: {row['pct_genes_k6plus']:.0f}% of genes have 6+ datasets. Disagreement is {i2_label(row['median_I2_k2plus']).lower()}, so check each gene's I².")
    if row["pct_genes_untested"] > 0:
        st.caption(f"{row['pct_genes_untested']:.0f}% of genes have no adjusted p-value (untested). Left out, not filled in.")

    # Volcano
    view = sub[sub["k_datasets"] >= min_k]
    if show == "Up in disease":
        view = view[view["direction"] == "up"]
    elif show == "Down in disease":
        view = view[view["direction"] == "down"]

    if view.empty:
        st.warning("No genes match these settings. Lower the minimum datasets or change the direction.")
    else:
        v = view.assign(
            direction=view["direction"].map({"up": "Up in disease", "down": "Down in disease"}),
            padj_text=view["padj_meta"].map(fmt_p),
            I2_text=view["I2"].map(lambda x: "n/a (one dataset)" if pd.isna(x) else f"{x:.0f}%"),
        )
        hover = {"gene_symbol": True, "log2FoldChange_meta": ":.2f", "padj_text": True, "I2_text": True, "k_datasets": True,
                 "neglog10_padj": False, "direction": False}
        common = dict(x="log2FoldChange_meta", y="neglog10_padj", symbol="direction", hover_data=hover,
                      symbol_map={"Up in disease": "triangle-up", "Down in disease": "triangle-down"},
                      labels={"log2FoldChange_meta": "log2 fold change (disease vs. comparison)",
                              "neglog10_padj": f"-log10 adjusted p-value (capped at {CAP})",
                              "padj_text": "adjusted p", "I2_text": "I²", "k_datasets": "datasets",
                              "gene_symbol": "gene"})
        title = f"{pretty(disease)}: the {len(v)} most significant genes"
        if multi:
            v = v.assign(log10_k=np.log10(v["k_datasets"]))
            fig = px.scatter(v, color="log10_k", color_continuous_scale="Viridis",
                             title=title + ". Color shows how many datasets support each gene (log scale)", **common)
            ticks = [t for t in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500) if v["k_datasets"].min() <= t <= v["k_datasets"].max()]
            fig.update_coloraxes(colorbar=dict(title="datasets<br>supporting", tickvals=np.log10(ticks), ticktext=[str(t) for t in ticks]))
        else:
            fig = px.scatter(v, color_discrete_sequence=[TIER_COLOR["low"]],
                             title=title + ". All from a single dataset, so color is not used", **common)
        fig.update_traces(marker=dict(size=9, opacity=0.8, line=dict(width=0.5, color="white")))
        fig.add_hline(y=-np.log10(0.05), line_dash="dot", line_color="#c0392b", annotation_text="padj = 0.05",
                      annotation_position="bottom right", annotation_font_color="#c0392b")
        top = float(v["neglog10_padj"].max())
        if n_at_cap := int(v["at_cap"].sum()):
            fig.add_hline(y=CAP, line_dash="dot", line_color="grey", annotation_text="display cap", annotation_position="top left")
        fig.update_yaxes(range=[0, (CAP if n_at_cap else top) * 1.1 + 1])
        fig.update_layout(height=520, title_font_size=15, legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, title=""),
                          margin=dict(t=90), plot_bgcolor="white")
        fig.update_xaxes(zeroline=True, zerolinecolor="lightgrey", gridcolor="#eeeeee")
        fig.update_yaxes(gridcolor="#eeeeee")
        st.plotly_chart(fig, use_container_width=True)
        n_cap = n_at_cap
        if n_cap:
            st.caption(f"{n_cap} genes sit on the display cap (padj reported as 0), so their order along the top means nothing.")

        with st.expander("Gene table and download"):
            table = v[["gene_symbol", "direction", "log2FoldChange_meta", "padj_meta", "I2", "k_datasets"]].rename(columns={
                "gene_symbol": "Gene", "direction": "Direction", "log2FoldChange_meta": "log2FC", "padj_meta": "Adjusted p",
                "I2": "I² (%)", "k_datasets": "Datasets"})
            st.dataframe(table, hide_index=True, use_container_width=True)
            st.download_button("Download these genes (CSV)", table.to_csv(index=False), f"{disease}_genes.csv", "text/csv")

    with st.expander("How to read this plot"):
        st.markdown("Right of zero = higher in disease. Higher = stronger evidence. **Dark = few datasets, yellow green = many.** "
                    "Shows the 150 strongest genes per direction.")

# ----------------------------------------------------------------------------- tab 3
with tab_compare:
    st.subheader("Where does a disease sit in the whole catalog?")
    pick = st.selectbox("Highlight a disease on the charts (optional)", ["(none)"] + options, format_func=lambda d: d if d == "(none)" else pretty(d))

    a, b = st.columns(2, gap="large")
    with a:
        bar = pd.DataFrame({"tier": TIER_ORDER, "diseases": tier_counts.values})
        bar["label"] = bar["tier"].map(TIER_LABEL)
        bar["share"] = bar["diseases"] / n_dis * 100
        fig_t = px.bar(bar, x="label", y="diseases", color="tier", color_discrete_map=TIER_COLOR,
                       text=bar.apply(lambda r: f"{r['diseases']} ({r['share']:.0f}%)", axis=1),
                       category_orders={"label": [TIER_LABEL[t] for t in TIER_ORDER]},
                       title=f"{pct_low:.0f}% of diseases have only a single dataset behind them",
                       labels={"label": "Confidence tier (best supported gene)", "diseases": "Number of diseases"})
        fig_t.update_traces(textposition="outside", marker_line_color="#333", marker_line_width=0.8)
        fig_t.update_layout(showlegend=False, height=420, title_font_size=15, plot_bgcolor="white", yaxis_range=[0, tier_counts.max() * 1.15])
        fig_t.update_yaxes(gridcolor="#eeeeee")
        if pick != "(none)":
            fig_t.add_annotation(x=TIER_LABEL[diseases.loc[diseases.disease == pick, "tier"].iloc[0]], y=0, yshift=-52,
                                 text=f"{pretty(pick)} is here", showarrow=False, font=dict(color="#c0392b", size=12))
        st.plotly_chart(fig_t, use_container_width=True)
    with b:
        rows = audit["rows"]
        g = pd.DataFrame({"tier": TIER_ORDER, "genes": [audit["k1_rows"], audit["k2_5_rows"], audit["k6_rows"]]})
        g["share"] = g["genes"] / rows * 100
        g["label"] = g["tier"].map(TIER_LABEL)
        fig_g = px.bar(g, x="share", y=["All gene results"] * 3, color="tier", color_discrete_map=TIER_COLOR, orientation="h",
                       text=g["share"].map(lambda x: f"{x:.0f}%"), hover_data={"genes": ":,", "label": True, "tier": False},
                       title=f"Counting every gene result, {audit['k1_rows'] / rows * 100:.0f}% rest on one dataset",
                       labels={"x": "Share of the 7.3 million gene results (%)", "share": "Share (%)", "y": ""})
        fig_g.update_traces(marker_line_color="#333", marker_line_width=0.8, textposition="inside",
                            textfont=dict(color="white"))
        fig_g.update_layout(barmode="stack", bargap=0.55, height=420, margin=dict(t=60), title_font_size=15, plot_bgcolor="white",
                            legend=dict(orientation="h", y=-0.25, title="", itemclick=False),
                            xaxis=dict(range=[0, 100], tickvals=[0, 25, 50, 75, 100], title="Share of all gene results (%)"), yaxis=dict(visible=False))
        for tr in fig_g.data:
            tr.name = TIER_LABEL[tr.name]
        st.plotly_chart(fig_g, use_container_width=True)
        st.caption("Left counts diseases, right counts gene results. Same story.")

    st.subheader("Do more datasets mean more agreement?")
    m = diseases[diseases["k_max"] >= 2].dropna(subset=["median_I2_k2plus"])
    rho, p_rho = stats.spearmanr(m["k_max"], m["median_I2_k2plus"])
    fig_s = px.scatter(m, x="k_max", y="median_I2_k2plus", color="tier", color_discrete_map=TIER_COLOR, size="n_sig", size_max=22,
                       log_x=True, hover_name=m["disease"].map(pretty),
                       hover_data={"k_max": True, "median_I2_k2plus": ":.0f", "n_sig": ":,", "tier": False},
                       title="Diseases with many datasets still show moderate to substantial disagreement between them",
                       labels={"k_max": "Datasets behind the best supported gene (log scale)",
                               "median_I2_k2plus": "Median I² of multi dataset genes (%)", "n_sig": "significant genes"})
    for y in (25, 50, 75):
        fig_s.add_hline(y=y, line_dash="dot", line_color="lightgrey")
    for tr in fig_s.data:
        tr.name = TIER_LABEL.get(tr.name, tr.name)
    if pick != "(none)" and pick in set(m["disease"]):
        r = m[m["disease"] == pick].iloc[0]
        fig_s.add_trace(go.Scatter(x=[r["k_max"]], y=[r["median_I2_k2plus"]], mode="markers+text", text=[pretty(pick)],
                                   textposition="top center", marker=dict(size=26, color="rgba(0,0,0,0)", line=dict(color="#c0392b", width=3)),
                                   showlegend=False, hoverinfo="skip"))
    fig_s.update_traces(marker_line_color="#333", marker_line_width=0.6, selector=dict(mode="markers"))
    fig_s.update_layout(height=480, title_font_size=15, plot_bgcolor="white",
                        legend=dict(orientation="h", y=1.02, x=0, yanchor="bottom", title=""), margin=dict(t=90))
    fig_s.update_xaxes(gridcolor="#eeeeee", tickvals=[2, 3, 5, 10, 20, 50, 100, 200, 500], ticktext=["2", "3", "5", "10", "20", "50", "100", "200", "500"],
                       range=[np.log10(1.6), np.log10(750)])
    fig_s.update_yaxes(gridcolor="#eeeeee", range=[0, 100], tickvals=[0, 25, 50, 75, 100],
                       ticktext=["0", "25 (low | moderate)", "50 (moderate | substantial)", "75 (substantial | considerable)", "100"])
    st.plotly_chart(fig_s, use_container_width=True)
    st.caption(f"{len(m)} diseases with 2+ datasets. Bubble size = significant genes. Spearman rho = {rho:.2f} (p = {p_rho:.3f}): weak upward trend.")
    if pick != "(none)" and pick not in set(m["disease"]):
        st.caption(f"{pretty(pick)} has only one dataset, so it cannot appear on this chart.")

    st.subheader("Do more datasets mean stronger signals?")
    gt = gene_sample.copy()
    gt["nlp"] = (-np.log10(gt["padj_meta"].clip(lower=1e-300))).clip(upper=CAP)
    rho_sig = stats.spearmanr(gt["k_datasets"], gt["nlp"])[0]
    rho_lfc = stats.spearmanr(gt["k_datasets"], gt["abs_lfc"])[0]
    med = gt.groupby("tier")["abs_lfc"].median().reindex(TIER_ORDER)
    fig_e = go.Figure()
    top_fence = 0
    for t in TIER_ORDER:
        g = gt.loc[gt["tier"] == t, "abs_lfc"]
        q1, q2, q3 = np.percentile(g, [25, 50, 75])
        iqr = q3 - q1
        hi, lo = g[g <= q3 + 1.5 * iqr].max(), g[g >= q1 - 1.5 * iqr].min()
        top_fence = max(top_fence, hi)
        fig_e.add_trace(go.Box(x=[TIER_LABEL[t]], q1=[q1], median=[q2], q3=[q3], lowerfence=[lo], upperfence=[hi], name=TIER_LABEL[t],
                               marker_color=TIER_COLOR[t], line=dict(color="#333"), fillcolor=TIER_COLOR[t], boxpoints=False, width=0.45))
    fig_e.update_layout(showlegend=False, height=380, title_font_size=15, plot_bgcolor="white",
                        title=f"Effect sizes shrink as datasets are added (median {med['low']:.2f} to {med['medium']:.2f} to {med['high']:.2f})",
                        yaxis_title="|log2 fold change|", xaxis_title="Datasets behind the gene")
    fig_e.update_yaxes(gridcolor="#eeeeee", range=[0, top_fence * 1.08])
    st.plotly_chart(fig_e, use_container_width=True)
    st.caption(f"{len(gt):,} random genes from all 350 diseases. Datasets vs significance: Spearman rho = {rho_sig:.2f} (no link). "
               f"Datasets vs |log2FC|: rho = {rho_lfc:.2f}. Box = middle 50%, line = median, whiskers = 1.5 x IQR. A larger first-study effect that shrinks when pooled is the usual 'winner's curse'.")

    st.subheader("Full disease list")
    f1, f2 = st.columns([1, 2])
    with f1:
        tiers_sel = st.multiselect("Confidence tier", TIER_ORDER, default=TIER_ORDER, format_func=lambda t: TIER_LABEL[t])
    with f2:
        text = st.text_input("Disease name contains", "", placeholder="for example: lung, carcinoma, diabetes")
    lst = diseases[diseases["tier"].isin(tiers_sel)]
    if text.strip():
        lst = lst[lst["disease"].str.replace("-", " ").str.contains(text.strip(), case=False, regex=False)]
    if lst.empty:
        st.warning("No diseases match. Clear the search box or select more tiers.")
    else:
        tbl = lst.sort_values("k_max", ascending=False)[["disease", "tier", "k_max", "k_median", "n_sig", "pct_genes_k6plus",
                                                          "median_I2_k2plus", "pct_genes_untested"]].copy()
        tbl["disease"] = tbl["disease"].map(pretty)
        tbl["Heterogeneity"] = tbl["median_I2_k2plus"].map(i2_label)
        tbl = tbl.rename(columns={"disease": "Disease", "tier": "Tier", "k_max": "Best gene datasets", "k_median": "Typical gene datasets",
                                  "n_sig": "Significant genes", "pct_genes_k6plus": "% genes with 6+ datasets",
                                  "median_I2_k2plus": "Median I² (%)", "pct_genes_untested": "% genes untested"})
        st.caption(f"{len(tbl)} of {n_dis} diseases. Click a column header to sort.")
        st.dataframe(tbl, hide_index=True, use_container_width=True,
                     column_config={"Significant genes": st.column_config.NumberColumn(format="%d"),
                                    "% genes with 6+ datasets": st.column_config.NumberColumn(format="%.1f"),
                                    "Median I² (%)": st.column_config.NumberColumn(format="%.0f"),
                                    "% genes untested": st.column_config.NumberColumn(format="%.1f"),
                                    "Typical gene datasets": st.column_config.NumberColumn(format="%.0f")})
        st.download_button("Download this list (CSV)", tbl.to_csv(index=False), "disease_confidence.csv", "text/csv")

# ----------------------------------------------------------------------------- tab 4
with tab_data:
    st.markdown(
        f"**Source:** OCTAD 2.0 meta analysis. **{n_dis} diseases, {audit['rows']:,} gene results.** "
        "One row = one gene in one disease."
    )
    st.subheader("Quality checks and decisions")
    na_k1 = audit["padj_na_k1"]
    checks = pd.DataFrame([
        ["Missing adjusted p", f"{audit['na']['padj_meta'] / audit['rows'] * 100:.1f}% of rows, all single dataset genes", "Not imputed. Gap is not random."],
        ["I² for one dataset", "Reported as 0 in every case", "Set to missing. Can't be measured."],
        ["padj exactly 0", f"{audit['padj_zero']:,} rows (raw p also 0)", f"Kept. Plotted capped at {CAP}."],
        ["|log2FC| above 10", f"{audit['abs_lfc_gt10']:,} rows, {audit['abs_lfc_gt10_k1'] / audit['abs_lfc_gt10'] * 100:.0f}% single dataset", "Kept. Color shows weak support."],
        ["Duplicate gene symbols", f"{audit['dup_symbol_rows'] / audit['rows'] * 100:.1f}% of rows", "Kept. No ID is clearly right."],
        ["Skewed values", "p-values near 0, counts from 1 to 577", "-log10 for p, log axes, 3 tiers."],
        ["Integration", "350 per disease files, same columns", "Stacked into one table. No keys to merge."],
        ["Encoding, scaling, imputation", "Disease is only a selector; nothing feeds a model", "N/A. Not needed for a descriptive app."],
    ], columns=["Check", "Found", "Decision"])
    st.table(checks.set_index("Check"))
    miss = pd.DataFrame({"Genes": ["1 dataset", "2+ datasets"], "pct": [na_k1 / audit["k1_rows"] * 100, 0.0]})
    fig_m = px.bar(miss, x="pct", y="Genes", orientation="h", text=miss["pct"].map(lambda x: f"{x:.1f}%"),
                   title="Missing p-values only occur in single dataset genes",
                   labels={"pct": "Share with no adjusted p-value (%)"}, color_discrete_sequence=[TIER_COLOR["low"]])
    fig_m.update_traces(textposition="outside", marker_line_color="#333", marker_line_width=0.8)
    fig_m.update_layout(height=220, title_font_size=15, plot_bgcolor="white", xaxis_range=[0, 20], margin=dict(t=50, b=20))
    st.plotly_chart(fig_m, use_container_width=True)
    with st.expander("Variables and data types"):
        st.table(pd.DataFrame([
            ["disease", "Categorical (350 levels)", "Selector and grouping"],
            ["gene_symbol", "Text", "Gene identity"],
            ["log2FoldChange_meta", "Continuous, already log scaled", "Effect size and direction"],
            ["padj_meta", "Continuous, 0 to 1, very skewed", "Significance (transformed)"],
            ["I2", "Continuous, 0 to 100, needs 2+ datasets", "Disagreement between datasets"],
            ["k_datasets", "Discrete count, 1 to 577", "Support (binned into 3 tiers)"],
        ], columns=["Variable", "Type", "Role"]).set_index("Variable"))
    with st.expander("Features I created"):
        st.markdown("- **Tier** from the best gene: 1, 2 to 5, 6+ datasets\n- **% significant of tested genes**\n"
                    "- **Median I²** of multi dataset genes\n- **-log10(padj)**, capped at 50")
    st.caption("Catalog counts use all 7.3M rows. Charts of distributions use a 3% random sample (seed 830).")

# ----------------------------------------------------------------------------- tab 5
with tab_limits:
    st.subheader("What this app cannot tell you")
    st.markdown(
        "- Dataset count is a proxy for confidence, not a measure of it.\n"
        "- `k_datasets` can reach the hundreds, so it may count comparisons, not studies.\n"
        "- Only the 300 strongest genes per disease are plotted.\n"
        "- Descriptive only. Tier cut points (1, 2 to 5, 6+) are my choice.\n"
        "- Tiers are partly confounded with disease type (well studied cancers vs rare conditions)."
    )
    st.markdown("**Other explanations to keep in mind**")
    st.markdown(
        "- Pooling more samples makes p-values smaller by construction, so high significance at high k is partly built in.\n"
        "- Shrinking effect sizes fit a winner's curse: a single small study tends to overstate an effect.\n"
        "- A low I² from two datasets is weak reassurance; two studies can barely show disagreement."
    )
    st.markdown("**Use it as a prompt to check the underlying datasets**, not to accept or discard a signature.")
    st.markdown("**Next step:** add a pathway view of the replicated genes, so a disease can be judged at the biology level, not only gene by gene.")
