# Mini Project: Benchmarking Boltz-2 Affinity Predictions on EGFR Kinase Inhibitors

**Question:** How well does Boltz-2 rank known EGFR inhibitors by potency, compared with simple cheminformatics baselines?

**Target:** Human EGFR (ChEMBL target `CHEMBL203`, UniProt `P00533`)

**Deliverable:** Reproducible Jupyter notebooks + README in a public GitHub repo.

## Environment

Windows 11, NVIDIA GeForce RTX 4060 Laptop GPU (8 GB VRAM), driver 581.86.
Miniconda, Python 3.11.

### Installation

I found during my installation that order matters. 
Installing boltz first can pull a CPU-only PyTorch build.

```bash
conda create -n boltz python=3.11 -y
conda activate boltz

# GPU-enabled PyTorch, from NVIDIA's index, not PyPI
pip install torch --index-url https://download.pytorch.org/whl/cu124
python -c "import torch; print(torch.cuda.is_available())"   # must print True

pip install "boltz[cuda]" -U
pip install rdkit pandas numpy scipy scikit-learn matplotlib seaborn \
            chembl_webresource_client py3Dmol pyyaml requests

# verify torch survived the boltz install
python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"
```

`environment.yml` is provided as a record of the exact resolved versions, but
`conda env create -f environment.yml` will install a CPU-only PyTorch, since the
custom index URL is not preserved. Follow the steps above instead. 
Use `conda env export --no-builds > environment.yml`

Python 3.13 does not work: boltz pins scipy==1.13.1, fails to build from source on Windows.

### Versions

- boltz: 2.0.3
- torch: 2.6.0+cu124
- rdkit: 2026.3.6

### Data Source

Raw Data: data\raw\CHEMBL203_activities.zip

Data source: ChEMBL 37, target CHEMBL203, downloaded 2026-09-24.
NB: The EGFR target report card: https://www.ebi.ac.uk/chembl/explore/target/CHEMBL203

### Processing raw input data

Use the notebook `notebook/chembl_data.ipynb` to process raw ChEMBL data to a clean .csv

Useful ChEMBL reference:
    - Activity data documentation: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions
    - webinar overview on ChEMBL database: https://www.ebi.ac.uk/training/events/chembl-workshop-drug-design/

Cleaning notes:
    - Removed invalid rows, e.g. missing pChEMBL data 
    - Only assays which measure binding affinity by IC50 (Ki, Kd are a small minority of binding affinity measurements)
    - Remove assays for mutants. Matches to strings in the "Assay Description"
        T790M 4732
        L858R 4773
        C797S 1399
        G719 1
        A431 1072
        L861Q 135
        E746 98
        \bexon 640
        \bdel 2354
        insert 97
        mutant 4523
        variant 0
    - Rows containing smiles strings which could not be standardized are removed. These failures are tracked in `data/processed/standardization_failures.csv`.
    - For duplicate measurements, the median pChEMBL Value is used in a cleaned data frame, saved at `data/processed/egfr_clean.csv`.
    - Median SD 0.5 log units across compounds with ≥3 measurements. This is essentially a noise floor for predictions by Boltz-2.

Summary of ChEMBL data filtering:
| Step                                                                                        |   N Rows |   Rows Removed |   N compounds |   Compounds Removed |
|:--------------------------------------------------------------------------------------------|---------:|---------------:|--------------:|--------------------:|
| Raw table                                                                                   |    58847 |            nan |         19118 |                 nan |
| Valid row : Binding Assay, Standard Relation, Homo Sapiens, pChEMBL Value, Valid            |    19866 |          38981 |         10761 |                8357 |
| Drop rows with missing smiles                                                               |    19866 |             -0 |         10761 |                  -0 |
| Limit the type of measurements to IC50, Ki and Kd                                           |    19544 |            322 |         10734 |                  27 |
| Assay Type, only IC50 for consistency of measure.                                           |    18123 |           1421 |         10424 |                 310 |
| Remove mutants, deletions and other variants.                                               |     9422 |           8701 |          6748 |                3676 |
| Standardized structures. See failed attempts at data/processed/standardization_failures.csv |     9414 |              8 |          6740 |                   8 |
| Merged duplicate measurements.                                                              |     6710 |           2704 |          6710 |                  30 |
| Removed non-organic entries.                                                                |     6709 |              1 |          6709 |                   1 |

### Exploratory Cheminformatics

Use the notebook `notebook/cheminformatics.ipynb` to explore cheminformatics on clean dataset. 

- Table of common descriptors + distribution plots + drug-like heuristics indicated, see `figures/molecular_descriptors.png`
- Lipinski RO5 compliance, only ~50% have 0 violations
- Mean ligand efficiency near 0.3, which is a desirable metric
- Most common Bemis–Murcko scaffolds account for ~15% of dataset, see `figures/top_scaffolds.png` -> diverse chemistry?
- Molecules classified as covalent, 4-Anilinoquinazoline-like, aminopurine-like, monocyclic or other, using SMARTS
- Whether using the Murcko scaffold or the manual classification, these 
groupings are strong predictors of potency. This is a lower bound of performance against which the model will be validated
| Classification            |   count |   median |     std |
|:--------------------------|--------:|---------:|--------:|
| 4-Anilinoquinazoline-like |    1577 |     7.31 | 1.17361 |
| Aminopurine-like          |     161 |     8.22 | 1.04396 |
| Covalent                  |    1838 |     7.07 | 1.27307 |
| Monocyclic                |      64 |     5.08 | 1.01039 |
| Other                     |    3069 |     6.26 | 1.25205 |
- Fingerprints were used to compute Tanimoto similarities on a random subset. Since the mean pairwise Tanimoto is only 0.16, the dataset is largely chemically diverse. See `figures\tanimoto_similarity_sample_heat.png` or `figures\tanimoto_similarity_sample_dist.png`
