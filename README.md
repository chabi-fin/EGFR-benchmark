# Mini Project: Benchmarking Boltz-2 Affinity Predictions on EGFR Kinase Inhibitors

**Question:** How well does Boltz-2 rank known EGFR inhibitors by potency, compared with simple cheminformatics baselines?

**Target:** Human EGFR (ChEMBL target `CHEMBL203`, UniProt `P00533`)

**Deliverable:** Reproducible Jupyter notebooks + README in a public GitHub repo.

**Key findings:**
1) On a structurally diverse set of 50 EGFR inhibitors, **Boltz-2's affinity
predictions ranked compounds less accurately** than simple cheminformatics
baselines. A nearest-neighbour lookup on Morgan fingerprints achieved Spearman
rho = 0.78 (95% CI [0.56,0.88]) against 0.23 (95% CI [-0.02,0.49]) for Boltz-2. The gap widened for novel chemistry: restricted to compounds with low similarity to the reference set, Boltz-2 showed no positive correlation (rho = -0.13, n = 17) while the fingerprint baseline retained rho = 0.49. *I suspect the poor performance of Boltz-2 on this target is related to the multiple binding modes of EGFR inhibitors.*
    
    > Note that the fingerprint baselines had access to measured potencies for
6,000+ EGFR compounds, while Boltz-2 predicted from structure alone. The
comparison reflects practical utility on a well-characterised target rather
than a like-for-like test of predictive capability.

2) **Boltz-2 potency predictions are compressed towards the mean**. Residuals correlate strongly with experimental potency (rho = -0.79, 95% CI [-0.87, -0.67]): the weakest compounds are overpredicted by a median of 1.1 log units and the most potent underpredicted by 1.7, with the model unbiased only around pChEMBL 6. This compression limits ranking performance directly.


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
python -m pip install torch --index-url https://download.pytorch.org/whl/cu124 --force-reinstall --no-cache-dir
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

- boltz: 2.2.1
- torch: 2.6.0+cu124
- rdkit: 2026.3.6

### Data Source

Raw Data: data\raw\CHEMBL203_activities.zip

Data source: ChEMBL 37, target CHEMBL203, downloaded 2026-09-24.
NB: The EGFR target report card: https://www.ebi.ac.uk/chembl/explore/target/CHEMBL203

## Notebooks

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

- Fingerprints were used to compute Tanimoto similarities on a random subset. Since the mean pairwise Tanimoto is only 0.16, the dataset is largely chemically diverse. See `figures/tanimoto_similarity_sample_heat.png` or `figures/tanimoto_similarity_sample_dist.png`
- In a principal component analysis over the fingerprints, the first two components explain ~11% of the variance. The classifiers cluster on PC1/PC2, see `figures\PCA_fingerprints_pred_class.png`. Additionally, the first two components very roughly separate the low affinity molecules from the rest of the dataset, see `figures/PCA_fingerprints_pred_activity.png`

### Subset selection

Use the notebook `notebooks/subset_selection.ipynb` to select a subset of molecules for activity predictions by Boltz-2. 

- Molecules are assigned to potency bins (low : 0-6, mid : 6-8, high : 8+)
- There is a moderate indication of potency, relative to manual classification.

| Potency   |   4-Anilinoquinazoline-like |   Aminopurine-like |   Covalent |   Monocyclic |   Other |
|:----------|----------------------------:|-------------------:|-----------:|-------------:|--------:|
| low       |                         206 |                  3 |        410 |           46 |    1269 |
| mid       |                         938 |                 64 |        991 |           18 |    1489 |
| high      |                         433 |                 94 |        437 |            0 |     311 |

- 15 molecules are selected per potency group. These selections are stratified by chemotype using rdkits's `MinMaxPicker` over the fingerprints (Morgan).
- 6 well-known drugs targeting EGFR are added to the subset. They span 4 binding types (type I, type II, type $1\tfrac{1}{2}$ and type VI (covalent); allosteric binding types are not included). A reference complex is available from the PDB for each binding type. Two of the drugs are covalent binders with no available experimental structure.
- Note: orininally, one additional covalent binder was included in the known references. Pelitinib was removed after matching the Clean InchiKey failed between the main and reference dataframes. This could indictate an issue with the protocol for standardization from Smiles. This will be ignored for now.
- Decoys/Inactives are added so the set includes approximate nonbinders. DUD-E decoys are sometimes identified as nonbinders using physicochemical properties alone. Instead, binders with *very* low activities are used as inactive decoys instead. These were matched on physicochemical properties. 
- Removed inactives (i.e. decoys) which are analogues of actives using fingerprint similarity, and selected 10 randomly for the subset
- Total of 61 molecules in the subset

### Prepare and run Boltz predictions

As above, use the notebook `notebooks/boltz_inputs.ipynb` to prepare and run inputs for Boltz-2.

- The human EGFR sequence was downloaded from UNIPROT (`P00533`). The sequence was truncated to the kinase domain [672–998], mimicking the numbering of the erlotinib complex (PDB: 1M17). Sequence alignment accounts for 1-indexing and signal peptide, and was verified against several key peptides.ref: Stamos et al. 2002 Structure.
- The clean smiles of each compound in the boltz subset and the kinase domain sequence were written to yaml files (total=61).
- A depiction of each molecule in the subset was drawn and summarized. See `figures/boltz_subset_mols/boltzs_subset_page[01-11].png
- Each compound + the kinase sequence were written to a yaml input file.
- The flags `--use_msa_server` and `--no_kernels` were used after troubleshooting a test prediction run. See notebook for details.
- Some predictions fail due to memory allocation errors. The batch run was restarted with the flag `--overwrite` until all the predictions were completed.
- The runtime for the subset predictions were written to `data/processed/runtime_log.csv`. This contains several outliers due to long idle time (CHEMBL3915508, t=76529 s) or restarted calculations on completed predictions (e.g. CHEMBL271410, 9.0s). Most computations needed ~5-15 minutes.
- Only one prediction computation failed, which was excluded from the subset (CHEMBL58)

### Evaluate Boltz predictions

Use the notebook `notebooks/evaluate_predictions.ipynb` to parse the Boltz-2 predictions.

- Affinity predictions and confidence estimations were extracted from Boltz prediction JSON files into a dataframe.
- The affinity predictions (IC50 [µM]) were converted to pIC50. 
- The table of compounds + prediction and confidence estimations are in `data/processed/boltz_preds.csv`
- The correlation between predicted and experimental pIC50 values were assessed using Pearson and Spearman correlations. Bootstrapping was used to generate a confidence interval for each, with Pearson: 0.34 with a confidence interval (CI) = [0.10, 0.60], p: 0.01 and Spearman: 0.23 CI = [-0.02, 0.49], p: 0.10. At this sample size (50), the data cannot distinguish a genuine moderate correlation to no correlation at all on compound ranking. 
- Scatter plot of prediciton vs. measurement (noncovalent). The uncertainty in the experimental values is shown in the shaded region. Color represents the Botlz-2 prediction on whether the ligand is a binder or nonbinder. 
![Scatter plot of prediciton vs. measurement (noncovalent)](figures/pred_vs_exp_covalent.png)
- ROC-AUC on the affinity probability binary with a cut off of 7 results in AUC 0.66 CI = [0.49, 0.80], (n=50, 16 positives). The lower bound sits at the random classifier values, suggesting weak discriminative ability. The result is not distinguishable from chance at this small sample size.   
![ROC-AUC for Boltz-2 binding prediction, $p\text{ChEMBL}_{\text{cut off}}=7$](figures/roc_auc.png)
- Boltz-2 beats the trivial benchmarks of heavy atoms (Spearman: 0.08) and molecular weight (Spearman: 0.06). Molecular size is a poor predictor of binding for this set, since the set was constructed by maximizing structural diversity and stratification into potency bins. 
- Boltz-2 is out performed by a QSAR baseline. A random forest on Morgan fingerprints, trained on the remaining ChEMBL data with the 40-compound subset held out, results in Spearman 0.76 CI = [0.59, 0.86]. Many compounds in the subset have close analogues in the training set. The Tanimoto similarity is only low within the subset. 
- Nearest Neighbor similarity outperforms Boltz-2 and matches performance to the random forest, Spearman 0.78 CI = [0.56, 0.88]. In this case, the pChEMBL of the test set is assigned the pChEMBL of the most Tanimoto-similar training compound. 
- Both Boltz-2 and Random forrest depend on chemical familiarity. Boltz-2 loses all its predictive power when there are no analogues in the training set. 
![Correlation depends on inclusion of analogues](figures/Boltz2_RF_Tanimoto_similarity.png)

### Conclusion: The structural modelling by Boltz-2 performs poorly compared to cheminformatics modelling based in substructures for the target EGFR. 

I would suspect that the structural modelling by Boltz-2 would give a predictive edge over cheminformatics / substructure modelling, *where the chemisty is novel*. However, that is not supported by the evaluation here. Any predictive power Boltz-2 may have for this target, disappears when the analogues are removed. Far cheaper cheminformatics baselines outperform Boltz-2, with and without analogues, for this target. A larger data set is needed for more conclusive statements. 

The low performance may be in part due to the diverse set of possible binding modes in EGFR. Others previously found Boltz-2 has low performance if the target lacks a well-defined binding pocket. **Therefore, I recommend against trusting Boltz-2 binding predictions in targets with multiple binding modes.** 
ref: [Wan et al. 2026 On the Reliability of AI Methods in Drug Discovery: Evaluation of Boltz-2 for Structure and Binding Affinity Prediction] 

While investigating whether Boltz-2 systematically underpredicts covalent bonds, I notices systematic regression towards the mean. There is clear negative correlation (Spearman -0.79, CI [-0.87,-0.67]) between experimental potency and the prediction residuals (predicted − experimental pChEMBL). **The weakest compounds are overpredicted by a median of one log unit, while the most potent are underpredicted by 1.7 log units**
![Correlation of experimental potency vs. residual potency](figures/exp_vs_residual_potency.png)
- The residuals decline monotonically across potency quadrents. The weakest compounds are overpredicted by a median of one log unit, while the most potent are underpredicted by 1.7.

| Experimental pChEMBL   |   n |   Median residual |
|:-----------------------|----:|------------------:|
| 4.0 – 5.3              |  18 |              1.06 |
| 5.3 – 6.7              |  16 |             -0.71 |
| 6.7 – 8.1              |  12 |             -0.83 |
| 8.1 – 9.4              |  14 |             -1.67 |
