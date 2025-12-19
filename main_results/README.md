# Main results: released embeddings and exploration scripts

This folder contains the **final ICD–ATC embeddings** learned on the ESND dataset
and a **minimal set of scripts** to explore their semantic structure.

The goal is to:
- share the embeddings used in the paper,
- provide simple, transparent examples for neighbourhood queries and latent-space visualisation.

No individual-level claims data are included.

---

## Contents

main_results/
├── README.md
├── data/
│ ├── embeddings_ESND_2FC.csv.gz # released embeddings
│ ├── vocab.csv # code ↔ index mapping
│ ├── label_table.csv # code ↔ label
│ ├── chapter_table.csv # code ↔ chapter (for colouring)
│ └── umap_chapter_subset.png # example UMAP projection
│
├── 01_neighbour_retrieval.R # nearest-neighbour queries
├── 02_umap_graph.R # UMAP and local graph visualisation
│
└── functions/
├── similarity_measures.R # cosine similarity
└── get_label.R # code → label helper



All scripts are intended to be run **from the repository root** and rely on **relative paths**.

---

## Data

### `embeddings_ESND_2FC.csv.gz`
Final static embeddings for ICD-10 diagnoses, ATC medications. Also contains french classiciations systems : procedures (CCAM), Devices (LPP) and Biology (NABM).

- Rows correspond to codes
- Columns correspond to embedding dimensions
- Code identities are provided via `vocab.csv`

### `vocab.csv`
Mapping between embedding rows and medical codes.

### `label_table.csv`
Human-readable labels for codes (used for interpretation).

### `chapter_table.csv`
High-level chapter assignment used for colouring UMAP projections.

---

## Scripts

### `01_neighbour_retrieval.R`
Simple helper to retrieve nearest neighbours of a target code
based on cosine similarity.

Typical use:
- inspect semantic neighbourhoods,
- explore diagnosis–treatment relations,
- sanity-check embedding coherence.

---

### `02_umap_graph.R`
Starter script for visual exploration of the embedding space.

Provides:
- a 2D UMAP projection of a selected subset of codes,
- colouring by chapter and shape by ontology,
- a small local similarity graph around a target code.

The UMAP parameters, subset selection and colour palette match those used in the paper.
The resulting figure is saved to `main_results/data/`.

---

## Notes

- The scripts are intentionally minimal and exploratory.
- They are not meant to be a general-purpose library.
- Users are expected to adapt code subsets, parameters and plots to their own needs.