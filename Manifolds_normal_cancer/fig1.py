# Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.ticker import AutoMinorLocator
from pathlib import Path
from scipy.stats.mstats import gmean
from sklearn.decomposition import PCA
from tqdm import tqdm

# Set output path and file names
output_path = Path('output')
sample_path = Path(r'../databases_external/TCGA/dataLIHC')
data_path = Path(r'../databases_external/TCGA/dataLIHC/data')
tissue = 'LIHC'

# Read sample data
sample = pd.read_excel(sample_path / 'sample.xls')

# Divide samples into normal and tumor groups
normal = sample['Sample Type'] == 'Solid Tissue Normal'
tumor = sample['Sample Type'] != 'Solid Tissue Normal'

# Read FPKM data and create a matrix of gene expression data
fpkm = pd.read_table(data_path / sample['File Name'][0], names=['gene', 'value'])
data = np.zeros([sample.shape[0], fpkm.shape[0]])
for i, filename in enumerate(tqdm(sample['File Name'], desc=f'Importing {tissue} data')):
    df = pd.read_table(data_path / filename, names=['gene', 'value'])
    data[i] = df['value']

# Normalize data and perform PCA
data = data + 0.1
ref = gmean(data[normal])
data = np.log2(data/ref)

U, S, VT = np.linalg.svd(data, full_matrices=0)
data_transform = data @ VT[:2].T

# Plot KDE plots for both normal and tumor groups
sns.kdeplot(x=-data_transform[normal, 0], y=-data_transform[normal, 1],
            cmap="Blues", fill=True)
sns.kdeplot(x=-data_transform[tumor, 0], y=-data_transform[tumor, 1],
            cmap="Reds", fill=True)

# Plot scatter plots for both normal and tumor groups
plt.scatter(-data_transform[normal, 0], -data_transform[normal, 1], c='b',
            edgecolors='white', linewidths=0.4, label='LIHC, normal')
plt.scatter(-data_transform[tumor, 0], -data_transform[tumor, 1], c='r',
            edgecolors='white', linewidths=0.35, marker='s', label='LIHC, tumor')
plt.legend(fontsize=13)

# Set plot limits, labels, and ticks
plt.xlabel('PC1', fontsize=15)
plt.ylabel('PC2', fontsize=15)
plt.xlim(-75, 300)
plt.ylim(-175, 150)

gca = plt.gca()
gca.xaxis.set_minor_locator(AutoMinorLocator(2))
gca.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.tick_params(axis='both', which='both', direction='in', labelbottom=True,
                labelleft=True, top=True, right=True, labelsize=13)
gca.xaxis.set_ticks_position('both')
gca.yaxis.set_ticks_position('both')

# Save and show plot
plt.tight_layout()
plt.savefig(output_path / 'fig1.pdf')
plt.savefig(output_path / 'fig1.png')
plt.show()
