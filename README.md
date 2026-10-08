# Models of site-specific breast cancer metastasis
Code and data used to train and evaluate models of distant and site-specific metastasis in “Molecular predictors of organ specific breast cancer metastasis” (Michmerhuizen et al., submitted for publication).

## Overview: 
Bulk tumor gene expression data and clinical site of first metastasis data from five publicly available breast cancer datasets (GSE2034, GSE2603, GSE12276, NKI295, METABRIC) were combined. Logistic regression models were trained using elastic net regularization to select gene expression features from primary tumors to classify patients based on metastatic potential and to further identify transcriptomic features contributing to site-specific metastasis. The scripts, data, and resulting models are included in this repository.

**[data](./data/)**: Processed clinical data and the expression datasets after batch correction using distance weighted discrimination (all datasets adjusted to METABRIC, which was used as reference). Expression data for each sample can be combined to generate the complete expression matrix which was too large for direct upload.

**[scripts](./scripts/)**: All scripts used to combine and process the original datasets and perform all analyses and generate figures. Sub-folders are numerically labeled in the order that they were run. 

**[models](./models/)**: Four models of metastasis highlighted in our manuscript: 
1. **Distant metastasis** (any site [brain, bone, liver, lung, multi-site])
    - 1a. Trained using gene expression signatures
    - 1b. Trained using individual genes
2. **Luminal bone metastasis**, trained using gene expression signatures
3. **Multi-site metastasis**, trained using gene expression signatures
