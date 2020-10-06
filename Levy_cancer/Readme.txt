The present directory is related to de work done 
in the paper arXiv:1507.08232, updated on July 2020

There are two files. A txt file named 'tissues_D-data.txt'
and a python script 'cancer.py'.

The script 'cancer.py' reads data on lifetime cancer risk,
number of normal stem cells (Nstem) and annual division rate per stem
cell (msc) for 31 different tissues. This data can be found on paper
https://doi.org/10.1126/science.1260825

Another set of data is imported, which includes the results obtained
by performing a principal component analysis (PCA) on gene expression
data downloaded from The Cancer Genome Atlas database for eight different
tissues. This includes centre of the tumour zone (x1) and the radii of the
normal zone (Rn) and tumour zone (Rt) along the first principal component (PC1).
This results can be found on paper
arXiv:1706.09813
and the calculations can be found on this repository
evolp/PCA_cancer/ 

A final set of data can be found on the txt file 'tissues_D-data.txt',
which collects the greatest modular contribution to the eigenvector 
associated to the first principal component (D) of the gene expression 
data for eight different tissues. The gene expression data was downloaded 
from The Cancer Genome Atlas (TCGA) database.
This analysis is performed on paper 
arXiv:2003.07828

Outputs for 'cancer.py'
- A data file with the values of x1, Rn, Rt, msc, Nstem, D, and cancer risk, see Table 1 on paper.
- A data file of the extra risk scores (ERS) for 31 different tissues, see Table 2 on paper.
- A plot of the downloaded samples for colon adenocarcinoma on the Gene Expression Space, see Fig.1 on paper.
- Two plots of ability to describe cancer risk of Eq.6 and Eq.10, see Fig.3 and Fig.4 on paper .
- A plot of cancer risk per stem cell vs time as total of stem cells divisions, see Fig.5 on paper.

Note: The cancer risk analysis is performed for the eight types of cancer with all parameters known.
COAD Colon adenocarcinoma
THCA Thyroid carcinoma
LUAD Lung adenocarcinoma
PRAD Prostate adenocarcinoma
BRCA Breast invasive carcinoma
LIHC Liver hepatocellular carcinoma
HNSC Head-Neck squamous cell carcinoma
ESCA Esophageal carcinoma
