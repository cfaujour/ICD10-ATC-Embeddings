# ICD-10–ATC Unified Embedding Space

This repository releases the **ICD-10–ATC Unified Embedding Space**:  
a set of joint, unsupervised embeddings for medical codes learned from large-scale healthcare claims data.

The **primary goal** of this repository is to **share these embeddings for reuse**
by researchers and data scientists working with EHR or claims data.

---

## Purpose and scope

The main interest of this resource is to provide a **shared latent space**
for **ICD-10 diagnosis codes** and **ATC medication codes**, enabling international researchers
to exploit semantic relationships between diagnoses and treatments beyond sparse, binary representations.

In addition, the embedding space also comprises **most medical coding systems used in the French SNDS**
(*Système National des Données de Santé*), including:
- **ICD-10** (diagnoses),
- **ATC** (medications),
- **CCAM** (procedures),
- **LPP** (medical devices),
- **NABM** (laboratory tests),

making it directly usable for French researchers working with SNDS-derived data.

No individual-level healthcare data are included.

---

## Intended usage

The embeddings are released as a **research resource**.
Typical uses include:

- **Automatic retrieval of related codes**  
  (e.g. semantic neighbourhoods, diagnosis–treatment associations, exploratory search).

- **Reuse or fine-tuning as an input embedding layer**  
  for machine-learning or neural models applied to claims or EHR data.

The scope is intentionally limited to these core use cases.

---

## Data source, method, and embedding properties

The embeddings are learned from longitudinal claims data extracted from the **SNDS**,
which covers nearly the entire French population and aggregates reimbursement and hospital information.

At a high level:
- Code sequences are extracted from care trajectories.
- Code–code co-occurrences are counted within temporal windows.
- The resulting matrix is transformed using **Positive Pointwise Mutual Information (PPMI)**.
- A **truncated Singular Value Decomposition (SVD)** is applied to obtain low-dimensional embeddings.

The resulting representations are:
- **unsupervised**,
- **static** (one vector per code),
- **task-agnostic**,
- shared across multiple medical coding systems,
- interpretable through neighbourhood analysis and low-dimensional projections.

---

## More information

Full methodological details, evaluation, and experimental results are described in the accompanying paper:

> *[Paper title]*  
> *NPJ Digital Medicine* (submitted)

---

## Repository structure

The repository is organised into two main subfolders.

### `main_results/` — Released embeddings and exploration starter kit

This is the **main entry point** for users interested in reusing the embeddings.

It contains:
- the released embeddings (CSV format),
- a `vocab.csv` file mapping **embedding row indices to medical code identifiers**
  (required to interpret or reuse the embeddings),
- auxiliary label and chapter tables for interpretation,
- lightweight scripts for:
  - nearest-neighbour queries,
  - latent-space visualisation (UMAP),
  - local similarity inspection.

Direct links to the released data:

- Embedding vectors:  
  `main_results/data/embeddings_ESND_2FC.csv.gz`

- Row index–code mapping:  
  `main_results/data/vocab.csv`

See [`main_results/README.md`](main_results/README.md) for details.

---

### `toy_data/` — Reproducible toy pipeline

This folder provides a **fully reproducible toy example** of the embedding pipeline,
using synthetic EHR-like data.

It illustrates each methodological step (co-occurrences, PPMI, SVD, visualisation)
for transparency and pedagogical purposes.

It is **not required** to reuse the released embeddings.

See [`toy_data/README.md`](toy_data/README.md).

---

## Software and usage notes

- The code is written in **R**.
- This repository is **not** a software package, API, or modelling framework.
- It is intended as a **data and method resource** rather than a production tool.

Cloning the repository is recommended to preserve the folder structure
and relative paths used in the scripts.

---

## Contact

For questions or reuse inquiries, please contact:  
[contact email]

---

## Citation

If you use the ICD-10–ATC Unified Embedding Space, or derived representations,
please cite the accompanying paper:

> *[Paper title]*  
> *Authors*  
> *NPJ Digital Medicine* (submitted)

This helps support continued sharing of research resources.

---

## Notes and potential extensions

This repository focuses on releasing embeddings and minimal exploration code.
Additional tools could be built on top of this resource.

In particular, a lightweight visualization interface allowing users to
interactively explore code neighbourhoods without writing code
(e.g. target-based semantic search and local projections)
would be a natural extension.

Contributions or suggestions in this direction are welcome.
