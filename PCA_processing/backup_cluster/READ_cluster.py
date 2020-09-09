import numpy as np
import xlrd as xl
from scipy.stats.mstats import gmean
from scipy.sparse.linalg import eigsh

def read_data():
    wb = xl.open_workbook("./../external_database/TCGA/20. TCGA-READ/sample.xls")
    ws = wb.sheet_by_index(0)

    directory = './../external_database/TCGA/20. TCGA-READ/'

    sample_type = []
    normal = []
    tumor = []
    data = []

    for i in range(1, ws.nrows):
        
        name = ws.cell_value(i,0)
        sample_type.append(ws.cell_value(i,3))

        file_object = open(file = directory + name, mode = "r")
        lines = file_object.readlines()

        if(ws.cell_value(i,3) == 'Solid Tissue Normal'):
            normal.append(i - 1)
        else:
            tumor.append(i - 1)

        data.append([])
        for line in lines:
            data[-1].append(line.split()[1])
        
        file_object.close()


    return([np.array(data, dtype='double') + 0.1, normal, tumor])

def calc(data, normal):
    ref = gmean(data[normal])

    t = np.log2(data/ref)

    covariance = np.dot(t.T, t)/t.shape[0]

    eigenvalues, eigenvectors = eigsh(covariance, k = 100)

    eigenvalues_normalized = eigenvalues/eigenvalues.sum()

    projection = np.dot(eigenvalues.T,t.T)
    
    return([eigenvalues, eigenvectors, eigenvalues_normalized, projection])

if __name__ == "__main__":
    data, normal, tumor = read_data()
    
    eigenvalues, eigenvectors, eigenvalues_normalized, projection = calc(data, normal)

    np.savetxt('./../PCA_tcga/READ/PC_dat.dat', projection.T, fmt='%f')
    np.savetxt('./../PCA_tcga/READ/eigenvectors_dat.dat', eigenvectors.T, fmt='%f')
    np.savetxt('./../PCA_tcga/READ/eigenvalues_normalized_dat.dat', eigenvalues_normalized, fmt='%f')
    np.savetxt('./../PCA_tcga/READ/eigenvalues_dat.dat', eigenvalues, fmt='%f')
    np.savetxt('./../PCA_tcga/READ/eigenvectorsT_dat.dat', eigenvectors, fmt='%f')

    principal = eigenvectors[:, eigenvalues.argmax()]
    index = np.argpartition(-np.abs(principal), 20)[:20]
    components = principal[index]

    np.savetxt('./../PCA_tcga/READ/20_index_dat.dat', index, fmt='%i')
    np.savetxt('./../PCA_tcga/READ/20_component_dat.dat', components, fmt='%f')
    np.savetxt('./../PCA_tcga/READ/normal_dat.dat', normal, fmt='%i')
    np.savetxt('./../PCA_tcga/READ/tumor_dat.dat', tumor, fmt='%i')