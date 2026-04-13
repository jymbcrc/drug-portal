from pathlib import Path
import streamlit as st

st.set_page_config(
    page_title="Project Overview",
    layout="wide"
)

BASE_DIR = Path(__file__).resolve().parent.parent

# =====================
# Title
# =====================
st.title("From proteomic screening to mechanistic insight with robotics and AI")

st.markdown("---")

# =====================
# Study overview
# =====================
st.header("Study overview")
st.markdown(
    '''
This interactive portal accompanies our study describing an end-to-end platform that integrates automated proteomic screening with AI-assisted data interpretation to enable scalable mechanism-of-action discovery.

The underlying dataset comprises a large-scale proteomic screen of 172 compounds in HepG2 cells, generating over 1,200 proteomes and quantifying more than 8,700 proteins. For each compound, protein-level differential analysis and pathway-level interpretation (Hallmark and Gene Ontology) were systematically performed, followed by AI-assisted summarization and hypothesis generation.

Rather than serving as a static data repository, this portal is designed to facilitate multi-level exploration of drug-induced cellular responses, enabling users to move seamlessly from global patterns to individual molecular details.

Key analytical layers available in this portal include:

- Protein-level regulation: differential abundance (log2FC, statistical significance) across all compounds  
- Pathway-level perturbation: GSEA-based enrichment across Hallmark and GO biological processes  
- Drug similarity structure: clustering and correlation-based grouping of compounds with shared responses  
- AI-generated molecular profiles: integrated summaries of drug function, regulated proteins, pathways, and potential novel insights  

Together, these features provide a structured framework for understanding both shared and compound-specific mechanisms of action.
'''
)

st.markdown("---")

# =====================
# How to use
# =====================
st.header("How to use this portal")
st.markdown(
    '''
This portal supports multiple complementary exploration modes:

1. Drug-centric exploration  
   Start from an individual compound to examine its proteomic signature, including differential proteins, volcano plots, and pathway perturbations.

2. Protein-centric exploration  
   Query a specific gene or protein to evaluate how it responds across the entire drug panel.

3. Structure discovery via clustering  
   - Identify groups of proteins with coordinated regulation patterns  
   - Identify clusters of drugs with similar cellular mechanisms  

4. Cross-scale interpretation  
   Move between protein-level, pathway-level, and drug-level views to build mechanistic insights.

This design enables both hypothesis-driven analysis and data-driven discovery.
'''
)

st.markdown("---")

# =====================
# Dataset summary
# =====================
st.header("Dataset summary")
st.markdown(
    '''
- 172 compounds screened in HepG2 cells  
- 1,232 proteomic samples across 13 plates  
- 6 biological replicates per compound  
- 8,703 quantified proteins  
- ~10 million data points  

Controls and quality design:
- DMSO (vehicle control)  
- Deferoxamine (DFO, positive control)  
- 12 QC samples per plate  

Analysis pipeline:
- Differential protein analysis (t-test with FDR correction)  
- Log2 fold-change–based ranking  
- Gene set enrichment analysis (Hallmark + GO BP)  
- AI-assisted interpretation and hypothesis generation  
'''
)

st.markdown("---")


# =====================
# Paper link
# =====================
st.header("Paper Link")
st.markdown(
"[Paper link](https://pmc.ncbi.nlm.nih.gov/articles/PMC10326973/)"
)

st.info(
    "This portal transforms large-scale proteomic data into an interactive, multi-layered resource for mechanism-of-action discovery and hypothesis generation."
)