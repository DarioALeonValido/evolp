The present directory is related to the work done 
in the paper arXiv:1507.08232, updated on July 2020.

There are two files. A txt file named 'tissues_D-data.txt'
and a python script 'cancer.py'.

The script 'cancer.py' reads data on lifetime cancer risk,
number of normal stem cells (Nstem) and annual division rate per stem
cell (msc) for 31 different tissues. This data can be found in the papers
https://doi.org/10.1126/science.1260825 and https://doi.org/10.1126/science.aaf9011.
It is also available in the databases_external/Cancer_Risk/ section of 
this repository.

A second set of data is imported, which includes the results obtained
by performing a principal component analysis (PCA) on gene expression
data downloaded from The Cancer Genome Atlas database for eight different
tissues. This includes the center of the tumor zone (x1) and the rms radii 
along the first principal component (PC1) of the normal zone (Rn) and 
tumor zone (Rt). These results can be found in the paper arXiv:1706.09813, 
updated in July 2020, and the way the calculations are performed may be 
found in this repository in the section PCA_cancer. 

A final set of data comes from the txt file 'tissues_D-data.txt',
which collects the greatest modular contribution (D) to the eigenvector 
associated to the first principal component PC1. See the paper
arXiv:2003.07828 for further details.

Outputs for 'cancer.py'
- A data file with the values of x1, Rn, Rt, msc, Nstem, D, and cancer risk, see Table 1 in the paper.
- A data file with the extra risk scores (ERS) for 31 different tissues, see Table 2 in the paper.
- A plot of the (PC1,PC2) plane for colon adenocarcinoma, see Fig.1 in the paper.
- Two plots of how Eq.6 and Eq.10 describe cancer risk of 8 tissues, see Fig.3 and Fig.4 in the paper .
- A plot of cancer risk per stem cell vs time, measured as number of stem cells divisions, see Fig.5 in the paper.

Note: The cancer risk analysis is performed for the following eight types of cancer with all parameters available.
COAD - Colon adenocarcinoma
THCA - Thyroid carcinoma
LUAD - Lung adenocarcinoma
PRAD - Prostate adenocarcinoma
BRCA - Breast invasive carcinoma
LIHC - Liver hepatocellular carcinoma
HNSC - Head and Neck squamous cell carcinoma
ESCA - Esophageal carcinoma
