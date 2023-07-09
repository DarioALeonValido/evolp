import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pathlib import Path
from scipy.stats.mstats import gmean
from sklearn.manifold import TSNE

sns.set(style="darkgrid")

# Set output path
output_path = Path('output')

# Set path to data
datapath = Path(r'..\databases_generated\TCGA_exp\meanvalues')

# Set tissue types
tissues = np.array(['BLCA', 'BRCA', 'COAD', 'ESCA', 'HNSC', 'KIRC', 'KIRP',
                    'LIHC', 'LUAD', 'LUSC', 'PRAD', 'READ', 'STAD', 'THCA',
                    'UCEC'])

# Load data
data = np.zeros((30, 60483))
for i, tissue in enumerate(tissues):
    data[i] = np.loadtxt(datapath / f'normal{tissue}.dat')
    data[i + 15] = np.loadtxt(datapath / f'tumor{tissue}.dat')

# Calculate reference (normal center)
ref = gmean(data[:15])

# Normalize data
data = np.log2(data/ref)

# Calculate tumor center
tumor_center = data[15:].mean(axis=0)
tumor_center_norm = tumor_center/np.linalg.norm(tumor_center)

# Calculate the distance between normal and its respective tumor
dist_tissues = data[15:] - data[:15]

# Create dataframe with results
df = pd.DataFrame(index=tissues)
df.index.name = 'Tissue'
df['Normalized projection'] = dist_tissues.dot(tumor_center_norm)/np.linalg.norm(dist_tissues, axis=1)
df['$\Theta(N)/\pi$'] = np.arccos(data[:15].dot(tumor_center_norm)/np.linalg.norm(data[:15], axis=1))/np.pi
dist = np.linalg.norm(data[15:] - data[:15, np.newaxis], axis=2)
df['Min(dist(N, T’))'] = dist.min(axis=1)
df['Closer Tumor'] = tissues[dist.argmin(axis=1)]
df['dist(N, T)'] = np.linalg.norm(data[15:] - data[:15], axis=1)

# Save results to files
open(output_path / 'table2.md', 'w').write(df.to_markdown())
