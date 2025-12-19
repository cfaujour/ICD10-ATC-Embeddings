# ICD-10–ATC Unified Embedding Space

This repository releases the **ICD-10–ATC Unified Embedding Space**, a set of joint, unsupervised vector representations (“embeddings”) for medical codes learned from large-scale healthcare claims data.  
Its primary purpose is to **share these embeddings for reuse** by researchers and data scientists working with EHR or claims data.

**Embeddings** represent each medical code as a dense numerical vector such that codes that tend to occur in similar clinical contexts (e.g. diagnoses and their associated treatments) are close to each other in the embedding space. This enables semantic comparisons between medical concepts beyond sparse, binary representations.

This resource is released in the context of the accompanying paper:

> *Medical Code Embeddings from Claims-Based Co-occurrences: A Unified Semantic Space for ICD-10 Diagnoses and ATC*  
> *NPJ Digital Medicine* (submitted) [link to be added]

---

## Quick access

- **Released embeddings:**  `main_results/data/embeddings_ESND_2FC.csv.gz`
- **Code vocabulary (row index–code mapping):**  `main_results/data/vocab.csv`

---

## Scope and code coverage

The embedding space provides a **shared latent representation**
for **ICD-10 diagnosis codes** and **ATC medication codes**, enabling international researchers
to explore semantic relationships between diagnoses and treatments within a unified vector space.

In addition, the embeddings also comprise **most medical coding systems used in the French SNDS**
(*Système National des Données de Santé*), including:
- **ICD-10** (diagnoses),
- **ATC** (medications),
- **CCAM** (procedures),
- **LPP** (medical devices),
- **NABM** (laboratory tests),

making the resource directly usable for French researchers working with SNDS-derived data.

---

## Intended usage

The embeddings are released as a **research resource**.
Typical uses include:

- **Automatic retrieval of related codes**  
  (e.g. semantic neighbourhoods, diagnosis–treatment associations, exploratory search).

- **Reuse as a pretrained embedding lookup table**
  optionally fine-tuned, for neural models applied to claims or EHR data using ICD-10, ATC, and related code systems.


---

## Data source, method, and embedding properties

The embeddings are learned from longitudinal claims data extracted from the **ESND** (subset of the SNDS),
which covers ~1.5 milion French care sequences and aggregates reimbursement and hospital information.

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

## Repository structure

The repository is organised into two main subfolders.

### `main_results/` — Released embeddings and exploration starter kit

This is the **main entry point** for users interested in reusing the embeddings.

It contains:
- the released embeddings,
- a `vocab.csv` file mapping **embedding row indices to medical code identifiers**
  (required to interpret or reuse the embeddings),
- auxiliary label and chapter tables for interpretation,
- lightweight scripts for:
  - nearest-neighbour queries,
  - latent-space (UMAP) and local graph visualisation

See [`main_results/README.md`](main_results/README.md) for details.

---

### `toy_data/` — Reproducible toy pipeline

This folder provides a **fully reproducible toy example** of the embedding pipeline,
using synthetic EHR-like data.

It illustrates each methodological step (co-occurrences, PPMI, SVD, visualisation)
for transparency and pedagogical purposes.

See [`toy_data/README.md`](toy_data/README.md).

---

## Software and usage notes

- The code is written in **R**.
- This repository is **not** a software package, API, or modelling framework.
- It is intended as a **data and method resource** rather than a production tool.

Cloning the repository is recommended to preserve the folder structure
and relative paths used in the scripts.

---

## Citation

If you use the ICD-10–ATC Unified Embedding Space, or derived representations,
please cite the accompanying paper:

> *Medical Code Embeddings from Claims-Based Co-occurrences: A Unified Semantic Space for ICD-10 Diagnoses and ATC *  
> *C.Faujour, S.Bouee, C.Emery, A-S.Jannot*  
> *NPJ Digital Medicine* (submitted) [link to come]

This helps support continued sharing of research resources.

---

## Contact

For questions or reuse inquiries, please contact:  
[corentin.faujour@cemka.fr]
[c.faujour@outlook.fr]

---

## Notes and potential extensions

This repository focuses on releasing embeddings and minimal exploration code.
Additional tools could be built on top of this resource.

In particular, a lightweight visualization interface allowing users to
interactively explore code neighbourhoods without writing code
(e.g. target-based semantic search and local projections)
would be a natural extension.

Contributions or suggestions in this direction are welcome.

