"""
MorphoChat: Morphological Learnability Atlas — Streamlit App

Which molecular signals are recoverable from H&E histology alone?
Assessed by spatially-aware block-permutation test with FDR correction.
"""
import json
from pathlib import Path
import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
    page_title="MorphoChat",
    page_icon="🔬",
    layout="wide",
)

@st.cache_data
def load_data(atlas_dir):
    atlases = {}
    for p in sorted(Path(atlas_dir).glob("learnability_*.json")):
        d = json.loads(p.read_text())
        atlases[d["organ"]] = d
    if not atlases:
        for p in sorted(Path(atlas_dir).glob("atlas_*.json")):
            if p.name == "atlas_summary.json":
                continue
            d = json.loads(p.read_text())
            atlases[d["organ"]] = d
    return atlases

atlas_dir = Path("atlas_data")
if not atlas_dir.exists():
    atlas_dir = Path("results/learnability")
if not atlas_dir.exists():
    st.error("Data not found.")
    st.stop()

atlases = load_data(str(atlas_dir))
has_perm = any("n_learnable_fdr" in a for a in atlases.values())
organs_sorted = sorted(atlases.keys(), key=lambda o: -(atlases[o].get("mean_observed_r") or atlases[o].get("summary", {}).get("mean_r") or 0))

total_genes = sum(a.get("n_genes", 200) for a in atlases.values())
total_learnable = sum(a.get("n_learnable_fdr", 0) for a in atlases.values())
total_slides = sum(a.get("n_slides", 0) for a in atlases.values())

st.title("🔬 MorphoChat")
st.markdown(f"""
**Which molecular signals are recoverable from H&E histology alone?**

Assessed by **spatially-aware block-permutation test** (100 permutations, BH-FDR < 0.05).
{total_slides} slides, {len(atlases)} organs, Virchow2 encoder.
**{total_learnable}/{total_genes}** gene-organ pairs have statistically significant morphological signal.

*Learnability is encoder-invariant (cross-encoder r > 0.95) — it's a property of the biology, not the model.*
""")

tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Explore Organ", "Search Gene", "About"])

with tab1:
    st.subheader("Learnability by Organ")
    rows = []
    for organ in organs_sorted:
        a = atlases[organ]
        if has_perm:
            nl = a.get("n_learnable_fdr", 0)
            ng = a.get("n_genes", 200)
            mr = a.get("mean_observed_r")
            mn = a.get("mean_null_r")
        else:
            s = a.get("summary", {})
            nl = s.get("n_learnable", 0)
            ng = nl + s.get("n_uncertain", 0) + s.get("n_not_learnable", 0)
            mr = s.get("mean_r")
            mn = None
        pct = nl / ng * 100 if ng > 0 else 0
        rows.append({
            "Organ": organ,
            "Slides": a.get("n_slides", 0),
            "Studies": a.get("n_studies", 0),
            "Spots": f'{a.get("n_spots", 0):,}',
            "Learnable (FDR)": f"{nl}/{ng} ({pct:.0f}%)",
            "Mean r": round(mr, 3) if mr else None,
            "Null r": round(mn, 4) if mn is not None else None,
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    st.subheader("Learnability Spectrum")
    bar_data = []
    for organ in organs_sorted:
        a = atlases[organ]
        if has_perm:
            nl = a.get("n_learnable_fdr", 0)
            ng = a.get("n_genes", 200)
            nn = ng - nl
        else:
            s = a.get("summary", {})
            nl = s.get("n_learnable", 0)
            nn = s.get("n_not_learnable", 0) + s.get("n_uncertain", 0)
        bar_data.append({"Organ": organ, "Learnable": nl, "Not significant": nn})
    bar_df = pd.DataFrame(bar_data).set_index("Organ")
    st.bar_chart(bar_df, color=["#2ecc71", "#e74c3c"])

with tab2:
    col1, col2 = st.columns(2)
    with col1:
        organ = st.selectbox("Organ", organs_sorted)
    with col2:
        if has_perm:
            verdict = st.selectbox("Filter", ["All", "Learnable (FDR<0.05)", "Not significant"])
        else:
            verdict = st.selectbox("Filter", ["All", "Learnable", "Uncertain", "Not Learnable"])

    a = atlases[organ]
    genes = a["genes"]

    if has_perm:
        if verdict == "Learnable (FDR<0.05)":
            genes = [g for g in genes if g.get("learnable", False)]
        elif verdict == "Not significant":
            genes = [g for g in genes if not g.get("learnable", False)]
    else:
        if verdict != "All":
            genes = [g for g in genes if g.get("verdict") == verdict.lower().replace(" ", "_")]

    gene_rows = []
    for g in genes:
        if has_perm:
            gene_rows.append({
                "Gene": g["gene"],
                "LOSO r": g.get("observed_r"),
                "Effect size (z)": g.get("effect_size"),
                "p-value": g.get("p_value"),
                "q-value (FDR)": g.get("q_value"),
                "Learnable": "Yes" if g.get("learnable") else "No",
            })
        else:
            gene_rows.append({
                "Gene": g["gene"],
                "LOSO r": g.get("r_loso"),
                "Self-conf": g.get("self_confidence"),
                "Stability": g.get("stability"),
                "Moran's I": g.get("morans_i"),
                "Verdict": g.get("verdict", "").replace("_", " ").title(),
            })

    st.markdown(f"**{organ}**: {a.get('n_slides', '?')} slides, {a.get('n_studies', '?')} studies, "
                f"{a.get('n_spots', 0):,} spots")
    st.dataframe(pd.DataFrame(gene_rows), use_container_width=True, hide_index=True, height=500)

with tab3:
    st.subheader("Is your gene learnable from H&E?")
    gene_query = st.text_input("Gene symbol", placeholder="e.g., MYL9, PODXL, COL1A2")

    if gene_query:
        gene_upper = gene_query.strip().upper()
        results = []
        for org, a in atlases.items():
            for g in a["genes"]:
                if g["gene"].upper() == gene_upper:
                    if has_perm:
                        results.append({
                            "Organ": org,
                            "LOSO r": g.get("observed_r"),
                            "Effect size (z)": g.get("effect_size"),
                            "q-value": g.get("q_value"),
                            "Learnable": "Yes" if g.get("learnable") else "No",
                        })
                    else:
                        results.append({
                            "Organ": org,
                            "LOSO r": g.get("r_loso"),
                            "Self-conf": g.get("self_confidence"),
                            "Verdict": g.get("verdict", "").replace("_", " ").title(),
                        })

        if results:
            st.dataframe(pd.DataFrame(results), use_container_width=True, hide_index=True)
            if has_perm:
                learnable_in = [r["Organ"] for r in results if r.get("Learnable") == "Yes"]
            else:
                learnable_in = [r["Organ"] for r in results if r.get("Verdict") == "Learnable"]
            if learnable_in:
                st.success(f"**{gene_upper}** is morphologically learnable in: {', '.join(learnable_in)}")
            else:
                st.warning(f"**{gene_upper}** is not significantly learnable in any tested organ")
        else:
            st.info(f"**{gene_upper}** not found in the top-200 variable genes for any organ.")

with tab4:
    st.subheader("Method")
    st.markdown("""
### Block-Permutation Learnability Test

For each gene in each organ, we test H₀: *"this gene's expression is not predictable from tissue morphology."*

**Procedure:**
1. Embed tissue patches (224×224 px) with frozen **Virchow2** (ViT-H) encoder
2. Select top-200 genes by within-organ expression variance
3. Compute LOSO Pearson r from Ridge regression on encoder features
4. Generate null distribution via **block permutation**: shuffle expression across spatial blocks (~25 spots) while preserving local autocorrelation
5. Compute empirical p-value from 100 permutations
6. Apply **Benjamini-Hochberg FDR correction** within each organ

**Key design choice:** Block permutation (not spot-level) preserves spatial autocorrelation under the null, preventing false positives from genes that appear predictable only because nearby spots have similar values.

**Result:** The null always produces r ≈ 0, confirming the block permutation destroys morphology-expression coupling. Effect sizes (z-scores) quantify learnability strength without arbitrary thresholds.

### Encoder Invariance
Cross-encoder per-gene learnability correlation (ResNet-50 vs UNI vs Virchow2) exceeds r = 0.79, reaching r = 1.00 between pathology FMs — confirming learnability is a biological property, not a model artifact.

**Data:** 111 HEST slides, 8 organs | **Encoder:** Virchow2 (Paige AI, ViT-H)
    """)
