import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from pathlib import Path
from scipy.stats.mstats import gmean

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

# Calculate angles between N and N'
n_vec = data[:15]/np.linalg.norm(data[:15], axis=1)[:, np.newaxis]
cos_nn = n_vec.dot(n_vec[:,:, np.newaxis])
angles_nn = np.arccos(np.squeeze(cos_nn)[np.eye(15, dtype=bool) == False])/np.pi

# Calculate angles between N and T'
t_vec = data[15:]/np.linalg.norm(data[15:], axis=1)[:, np.newaxis]
cos_nt = n_vec.dot(t_vec[:,:, np.newaxis])
angles_nt = np.arccos(cos_nt.flatten())/np.pi

# Calculate angles between T and T' and difference with angles between N and N'
data_t = data - tumor_center
tt_vec = data_t[15:]/np.linalg.norm(data_t[15:], axis=1)[:, np.newaxis]
cos_tt = tt_vec.dot(tt_vec[:,:, np.newaxis])
angles_tt = np.arccos(np.squeeze(cos_tt)[np.eye(15, dtype=bool) == False])/np.pi
diff = angles_tt - angles_nn

# Create a grid of panels to plot the distributions
fig = plt.figure(constrained_layout=True, figsize=(10, 10))
gs = fig.add_gridspec(2, 4)
ax1 = fig.add_subplot(gs[0,:2])
ax2 = fig.add_subplot(gs[0,2:])
ax3 = fig.add_subplot(gs[1,1:3])

fontsize_axes = 14
fontsize_title = 16
fontsize_ticks = 13
plt.tick_params(labelsize=fontsize_ticks)

# Plot the distribution of angles_nn in the first panel
sns.histplot(angles_nn, kde=False, stat="frequency", ax=ax1)
sns.rugplot(angles_nn, ax=ax1, height=0.06)
ax1.set_ylabel('Frequency', fontsize=fontsize_axes)
ax1.set_xlabel('Angle/Pi', fontsize=fontsize_axes)
ax1.set_title('a) Angle between N and N\'', fontsize=fontsize_title)

# Plot the distribution of angles_nt in the second panel
sns.histplot(angles_nt, kde=False, stat="frequency", ax=ax2)
sns.rugplot(angles_nt, ax=ax2, height=0.06)
ax2.set_ylabel('Frequency', fontsize=fontsize_axes)
ax2.set_xlabel('Angle/Pi', fontsize=fontsize_axes)
ax2.set_title('b) Angle between N and T\'', fontsize=fontsize_title)

# Plot the distribution of diff in the third panel
sns.histplot(diff, kde=False, stat="frequency", ax=ax3)
sns.rugplot(diff, ax=ax3, height=0.06)
ax3.set_ylabel('Frequency', fontsize=fontsize_axes)
ax3.set_xlabel('Angle/Pi', fontsize=fontsize_axes)
ax3.set_title('c) Difference between angles T-T\' and N-N\'', fontsize=fontsize_title)

# Save the figure to a file
fig.savefig(output_path / 'figure_3.png')
fig.savefig(output_path / 'figure_3.pdf')
plt.show()
