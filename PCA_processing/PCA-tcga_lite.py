import numpy as np
import xlrd as xl
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
from pathlib import Path
import PCA_core as pca

# Names ----------------------------------------------------------

lite_version = True
Ngen = 1000
tissue_id = "GBM"
#possible targets:
#tissue_id=["BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC",
#                  "PRAD","STAD","THCA","UCEC","BLCA","ESCA","READ"]

datapath = "../external_databases/TCGA/"
outputpath = "../PCA_tcga/"

# Functions ------------------------------------------------------

def read_data(datapath,tissue_id):
    wb = xl.open_workbook(datapath+"samples"+tissue_id+".xls")
    ws = wb.sheet_by_index(0)
    print(ws.nrows,"individual files were identified.")

    sample_type = []
    normal = []
    tumor = []
    data = []

    #reading first file to evaluate the cost
    name = ws.cell_value(1,0)
    sample_type.append(ws.cell_value(1,3))
    file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
    lines = file_object.readlines()
    print("Each file contains",len(lines),"genes.")
    print("Loading large databases...")
    if(ws.cell_value(1,3) == 'Solid Tissue Normal'):
        normal.append(0)
    else:
        tumor.append(0)
    data.append([])
    for line in lines:
        data[-1].append(line.split()[1])
    file_object.close()

    #reading the rest of the files
    for i in range(2, ws.nrows):
        name = ws.cell_value(i,0)
        sample_type.append(ws.cell_value(i,3))
        file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
        lines = file_object.readlines()
        if(ws.cell_value(i,3) == 'Solid Tissue Normal'):
            normal.append(i - 1)
        else:
            tumor.append(i - 1)
        data.append([])
        for line in lines:
            data[-1].append(line.split()[1])
        file_object.close()

    data = np.array(data, dtype="f") + 0.1
    ref = gmean(data[normal])  ## esto no me funciono asi
    # data_n = []
    # for i in normal:
    #     data_n.append(data[i])
    # ref = gmean(data_n)

    data = np.log2(data/ref)

    return([data, normal, tumor])



# Reading and processing the data ------------------------------------------

print("Reading files from databases:")
#if __name__ == "__main__":
data, normal, tumor = read_data(datapath,tissue_id)
print("Data successfully loaded!")
print("normal =", len(normal))
print("tumor =", len(tumor))
print("total =", len(data))

# reducing the data by 30 genes for lite version
if(lite_version and data.shape[1] > Ngen):
    print("Lite version: reducing the number of genes by", Ngen)
    data = data[:, :Ngen]
elif(len(data) < Ngen):
    lite_version = False
    print("The number of genes selected for lite version can't be upper that the number of genes in the files")
    print("Returning to the full version")
# data_red = []
# for i in range(len(data)):
#     data_red.append([])
#     for j in range(30):
#         data_red[i].append(data[i][j])

print("Computing PCA components...")
#eigenvalues, eigenvectors, eigenvalues_normalized, projection = pca.calc(data, normal) #Full version
eigenvalues, eigenvectors, eigenvalues_normalized, projection = pca.PC_decomposition(data, normal)
principal = eigenvectors[:, eigenvalues.argmax()]
index = np.argpartition(-np.abs(principal), 20)[:20]
components = principal[index]
print("Done!")


# Output data and figures --------------------------------------------------

plt.grid(True)
plt.xlabel('PC1')
plt.ylabel('PC2')

plt.scatter(projection[0,normal], projection[1,normal], c='b', s=15,label="Normal")
plt.scatter(projection[0,tumor], projection[1,tumor], c='r', s=15,label="Tumor")

plt.legend()


print("Exporting data...")
outputpath = outputpath + tissue_id + '_results/'
Path(outputpath).mkdir(exist_ok=True)
np.savetxt(outputpath + 'PC_dat.dat', projection, fmt='%f')
np.savetxt(outputpath + 'eigenvectors_dat.dat', eigenvectors.T, fmt='%f')
np.savetxt(outputpath + 'eigenvalues_normalized_dat.dat', eigenvalues_normalized, fmt='%f')
np.savetxt(outputpath + 'eigenvalues_dat.dat', eigenvalues, fmt='%f')
np.savetxt(outputpath + 'eigenvectorsT_dat.dat', eigenvectors, fmt='%f')
np.savetxt(outputpath + '20_index_dat.dat', index, fmt='%i')
np.savetxt(outputpath + '20_component_dat.dat', components, fmt='%f')
np.savetxt(outputpath + 'normal_dat.dat', normal, fmt='%i')
np.savetxt(outputpath + 'tumor_dat.dat', tumor, fmt='%i')

plt.savefig(outputpath + tissue_id + '_PC2_Vs_PC1_fig.png')
plt.show()

print("Done!")
