# Signature Confidence Explorer

An interactive Streamlit app that helps a researcher judge **how far to trust an OCTAD 2.0 disease gene signature**. For each of the 350 diseases it shows how many independent datasets support the signature and how much those datasets agree.

**Main finding:** 61% of the 350 diseases (213) rest on a single dataset, so their signatures cannot be checked for replication. Only 35% (123) have a gene backed by 6 or more datasets, and diseases with more than one dataset show moderate to substantial disagreement between them (median I² 53%).

> Live app: https://cmse830-signature-confidence.streamlit.app/

## Who it is for

A researcher about to build on a disease signature (for example, comparing diseases or choosing a disease model) who asks: *is this signature backed by more than one study, and do those studies agree?* The app describes the evidence. It does not recommend drugs or test hypotheses.

## What is in the app

| Tab | Purpose |
|---|---|
| Start here | Who it is for, the main message, how to use it, plain language definitions |
| Explore a disease | Confidence card (tier, datasets, significant genes, I² label) and a volcano plot of the strongest genes colored by dataset support, with filters and CSV download |
| Compare diseases | Catalog wide imbalance, dataset count vs I², and a filterable, sortable list of all 350 diseases |
| Data and preparation | Source, observation unit, every quality check with counts and the decision made |
| Limitations | What the app cannot tell you, other explanations, next step |

## Where to find each rubric area

| Rubric area | Where |
|---|---|
| Problem, audience, narrative | App tab *Start here*; notebook 02 section 1 |
| Data quality and preparation | App tab *Data and preparation*; notebook 02 sections 2, 3 and 7; notebook 01 |
| EDA and visualization | App tabs *Explore a disease* and *Compare diseases*; notebook 02 sections 4 to 6b |
| Interpretation and limitations | Notebook 02 interpretation blocks and section 8; app tab *Limitations* |
| Streamlit data product | `app.py` (5 tabs, 6 interactive controls) |
| GitHub and reproducibility | This README, `requirements.txt`, notebooks 01 (rebuild the data) and 02 (reproduce the figures), `data/processed/` |
| AI disclosure | Below; also in notebook 02 |

## Run it

```bash
git clone https://github.com/shubham-koirala/CMSE830.git
cd CMSE830/Midterm/signature-confidence-explorer
pip install -r requirements.txt
streamlit run app.py
```

The processed data (about 21 MB) is included in data/processed/, so the app runs as is. To open the notebooks, also pip install jupyter matplotlib. Notebook 01 needs the raw file, which is not shared, but its saved outputs show every step.

## Data and data access

- **Source:** OCTAD 2.0 random effects meta analysis signatures, one `<Disease>_meta.csv` per disease (350 files, about 7.3 million gene rows). Columns used: `log2FoldChange_meta`, `pval_meta`, `padj_meta`, `I2`, `k_datasets`, `gene_symbol`, `se_meta`.
- **Included in this repo:** the processed tables in `data/processed/` (about 21 MB). They are all the app needs, so the app runs without the raw data.
- **Not included:** the raw data (1.8 GB combined) is too large to share here. Notebook `01_prepare_data.ipynb` documents exactly how the processed tables were built from it (full pass over all 7.3 million rows, 3% random sample with seed 830, top 150 genes per direction per disease), and its outputs are saved in the notebook.

| File in `data/processed/` | Content |
|---|---|
| `disease_summary.csv` | one row per disease: datasets, significant genes, I², % untested, tier |
| `top_genes.csv` | 150 strongest up and 150 strongest down genes per disease |
| `gene_sample.csv` | 3% random sample of all gene rows (seed 830), used for catalog wide statistics |
| `audit.json` | exact counts from all 7.3 million rows, used in the Data tab and notebook |

## Repository layout

```
app.py                      Streamlit app
notebooks/01_prepare_data.ipynb   raw export -> small tables
notebooks/02_eda_and_quality.ipynb data quality audit, decisions and EDA (executed)
data/processed/             app ready data
docs/Signature_Confidence_Explorer_Midterm.pdf
requirements.txt, .streamlit/config.toml
```

## Key decisions (details in the app and notebook)

- Tier = datasets behind a disease's best gene: low (1), medium (2 to 5), high (6 or more).
- Missing adjusted p-values (8.1% of rows, all single dataset genes) are not imputed.
- I² is set to missing for single dataset genes, where it cannot be measured.
- -log10 adjusted p-value is capped at 50 for display, genes on the cap are not ranked against each other.
- No scaling, encoding or modeling: the project is descriptive. Modeling is out of scope for this project.

## Limitations

Dataset count is a proxy for confidence, not a measurement of it, and `k_datasets` can reach the hundreds, so it may count comparisons rather than distinct studies. The app plots only the strongest 300 genes per disease. Tier cut points are my choice. Differences between tiers are partly confounded with disease type.

## AI use disclosure

- **Tool and version:** Claude Sonnet 5.5
- **Purposes:** Drafting python scripts for plotting code, the Streamlit app
- **Representative prompts:** Asked few things during the ideation process such as: "What are the pros and cons of using the whole OCTAD dataset?"
- **How outputs were checked:** All counts were recomputed from the full 7.3 million row file and compared with the app and this notebook, the earlier hand picked sample for HW5 was replaced after I noticed it gave a misleading tier split, chart scales, transformations and labels were inspected.
- **My own contribution:** the project question, framing, decision, the choice to treat dataset support as the confidence signal, the decisions in section 3, and interpretations.