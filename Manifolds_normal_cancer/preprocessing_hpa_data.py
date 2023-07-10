import numpy as np
import pandas as pd

from functools import reduce
from pathlib import Path

datapath = Path("../databases_external/HPA/")
output_path = Path('../databases_generated/HPA_exp/')
v = np.array([1, 0, -1, -1])
values = {  # equivalence for data values
    'High': v[0],
    'Medium': v[1],
    'Low': v[2],
    'Not detected': v[3]
}


# Read data files
tumor = pd.read_table(datapath / 'pathology.tsv.zip')
normal = pd.read_table(datapath / 'normal_tissue.tsv.zip')

# Select only tissues with over 10,000 genes
tissues_to_keep = normal.groupby('Tissue').Gene.nunique()
tissues_to_keep = tissues_to_keep[tissues_to_keep >= 10000].index

normal = normal[normal.Tissue.isin(tissues_to_keep)]

# Reset the index
normal = normal.reset_index(drop=True)

# Select only the rows that have values for all four columns
tumor = tumor[tumor[['High', 'Medium', 'Low', 'Not detected']].notna().all(axis=1)]

# Eliminate genes that are not associated with all remaining cancer types
nc = tumor.Cancer.nunique()
gene_to_keep = tumor.groupby('Gene').Cancer.nunique()
gene_to_keep = gene_to_keep[gene_to_keep == nc].index

tumor = tumor[tumor.Gene.isin(gene_to_keep)]

# Reset the index
tumor = tumor.reset_index(drop=True)

# Obtaining all the tissues and types of cancer
cancer = tumor.Cancer.unique()
tissue = normal.Tissue.unique()

normal_genes = normal.groupby('Tissue').Gene.unique()
tumor_genes = tumor.groupby('Cancer').Gene.unique()
genes = reduce(np.intersect1d, pd.concat([normal_genes, tumor_genes]))

normal = normal[normal.Gene.isin(genes)]
normal = normal[normal.Level.isin(values.keys())]
normal = normal.replace({'Level': values})
normal.Level = normal.Level.astype(float)
normal = normal.groupby(['Tissue', 'Gene']).Level.mean()
normal = normal.unstack(level=-1).to_numpy()


tumor = tumor[tumor.Gene.isin(genes)]
total = tumor.groupby(['Cancer', 'Gene'])[['High', 'Medium', 'Low', 'Not detected']].sum().sum(axis=1)
tumor.loc[:, 'High'] = tumor['High'].astype(float) * v[0]
tumor.loc[:, 'Medium'] = tumor['Medium'].astype(float) * v[1]
tumor.loc[:, 'Low'] = tumor['Low'].astype(float) * v[2]
tumor.loc[:, 'Not detected'] = tumor['Not detected'].astype(float) * v[3]
tumor = tumor.groupby(['Cancer', 'Gene'])[['High', 'Medium', 'Low', 'Not detected']].sum().sum(axis=1) / total
tumor = tumor.unstack('Gene').to_numpy()

# Save data
np.savetxt(output_path / 'normalHPA.dat', normal)
np.savetxt(output_path / 'tumorHPA.dat', tumor)
