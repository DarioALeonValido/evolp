"""
This script file computes the entropy and overlap of PCA. 
The procedure and main results are described in paper 
arXiv:2005.02271
For more details on this file see author(s):
FEQ
"""

import numpy as np
import math, pylab
import xlrd as xl
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
import sys,os


# General variables ------------------------------------------------------

srcpc = "../databases_generated/TCGA_pca/"
srcsampl = "../databases_external/TCGA/"

# Functions ------------------------------------------------------
def read_data(srcpc,srcsampl):
    data={}
    pcfiles = os.listdir(srcpc)
    samplefiles = os.listdir(srcsampl)
    samplname = [f[+6:-4] for f in samplefiles if f[-3:] == 'xls']
    pcname = [f[+2:-4] for f in pcfiles if f[-3:] == 'xls']
    missmatch_pc_samples=list(set(pcname) - set(samplname))
    for name in samplname:
        if (name not in missmatch_pc_samples):
            data[name]={}
            data[name]['normal'], data[name]['tumor']=[],[]
            wbpc,wbsampl = xl.open_workbook(srcpc+"pc"+name+".xls"),xl.open_workbook(srcsampl+"sample"+name+".xls")
            wspc,wssampl = wbpc.sheet_by_index(0),wbsampl.sheet_by_index(0)
            minrows=min(wssampl.nrows,wspc.nrows)
            for row in range(1,minrows):
                if (wssampl.cell_value(row,3)=='Solid Tissue Normal'):
                    data[name]['normal'].append(wspc.row_values(row-1))
                else:
                    data[name]['tumor'].append(wspc.row_values(row-1))
    return data, missmatch_pc_samples

def _transpose_data(data):
    core_matrix = {}
    for label in data.keys():
        core_matrix[label] = {}
        for stage, core_matrix_T in data[label].items():
            core_matrix[label][stage] = np.around(np.array(core_matrix_T, dtype=np.float64).T, decimals = 10)
    return core_matrix

def _low_dimension_core_matrix(core_matrix, LOW_DIMENSION):
        unif_restricted_core_matrix={}
        for label in core_matrix.keys():
            unif_restricted_core_matrix[label]={}
            for stage in core_matrix[label].keys():
                unif_restricted_core_matrix[label][stage]=[]
                for i in range(LOW_DIMENSION):
                    row=core_matrix[label][stage][i]
                    unif_restricted_core_matrix[label][stage].append(row)
        return unif_restricted_core_matrix

def _PC_pair_plot(core_matrix, i, j,location_name):
    PC_j={}
    PC_i={}
    names = []
    for label in core_matrix.keys():
        names.append(label)
        PC_i[label]={}
        PC_j[label]={}
        for stage in core_matrix[label].keys():
            PC_i[label][stage]=core_matrix[label][stage][i-1]
            PC_j[label][stage]=core_matrix[label][stage][j-1]
        
    fig, ax = plt.subplots()
    ax.scatter(PC_i[location_name]['normal'], PC_j[location_name]['normal'], label='normal')
    ax.scatter(PC_i[location_name]['tumor'],PC_j[location_name]['tumor'], label='tumor')
    
    plt.xlabel("pc_{0}".format(i))
    plt.ylabel("pc_{0}".format(j))
    plt.title("pc_map_{0}".format(location_name))
    pylab.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig('pairplot_fig.pdf')   
    return None

def _std_dict(core_matrix):
    std={}
    for label in core_matrix.keys():
        std[label]={}
        for stage in core_matrix[label].keys():
            std[label][stage]=[]
            for row in core_matrix[label][stage]:
                std[label][stage].append(np.std(row))
    return std

def _mean_dict(core_matrix):
    mean={}
    for label in core_matrix.keys():
        mean[label]={}
        for stage in core_matrix[label].keys():
            mean[label][stage]=[]
            for row in core_matrix[label][stage]:
                mean[label][stage].append(np.mean(row))
    return mean

def _covariance_matrix_(core_matrix):
    covariance_matrix={}
    for label in core_matrix.keys():
        covariance_matrix[label]={}
        for stage in core_matrix[label].keys():
            covariance_matrix[label][stage]=\
            np.cov(core_matrix[label][stage], ddof=1 )
    
    return covariance_matrix

def _matrix_SVD(covariance_matrix):
    s, cov_inv = {},{}
    for label in covariance_matrix.keys():
        s[label], cov_inv[label] = {},{}
        for stage in covariance_matrix[label].keys():
            u, ss, vh =np.linalg.svd(covariance_matrix[label][stage])
            s[label][stage]=np.diag(ss)
            cov_inv[label][stage]=np.matmul(np.matmul(vh.T,np.linalg.inv(s[label][stage])),u.T)
    return s, cov_inv

def _H_multiv_gaussian(covariance_matrix,DIMENSION):
    H_={}
    S_, cov_inv = _matrix_SVD(covariance_matrix)
    for label in covariance_matrix.keys():
        H_[label]={}
        for stage in covariance_matrix[label].keys():
            sum_log = sum(math.log(S_[label][stage][i][i]) for i in range(DIMENSION))
            H_[label][stage] =\
            0.5*(sum_log + DIMENSION*(1 + math.log(2*math.pi) ))   
    return H_

def _delta_H(entropy_dict):
    delta_H={}
    for label in entropy_dict.keys():
        delta_H[label]= entropy_dict[label]['tumor']\
                       -entropy_dict[label]['normal']
    return delta_H                   

def _overlap(mean, covariance_matrix, DIMENSION):
    overlap, log_overlap={},{}
    S, Cov_Inv=_matrix_SVD(covariance_matrix)
    A, B = DIMENSION/2 , DIMENSION/4
    for label in mean.keys():
        inv_Vc = Cov_Inv[label]['normal']+Cov_Inv[label]['tumor']
        Vc = np.linalg.inv(inv_Vc)
        
        u_Vc, s_Vc, vh_Vc = np.linalg.svd(Vc)
        u_inv_Vc, s_inv_Vc, vh_inv_Vc = np.linalg.svd(inv_Vc)
        det_Vc, det_inv_Vc = np.prod(s_Vc), np.prod(s_inv_Vc)
        
        u_inv_Vn, s_inv_Vn, vh_inv_Vn = np.linalg.svd(Cov_Inv[label]['normal'])
        u_inv_Ve, s_inv_Ve, vh_inv_Ve = np.linalg.svd(Cov_Inv[label]['tumor'])
        
        
        log_sum = sum( -math.log(s_inv_Vn[i]) - math.log(s_inv_Ve[i]) +math.log(s_inv_Vc[i])  for i in range(DIMENSION))
        
        
        eta_n = np.matmul(Cov_Inv[label]['normal'],mean[label]['normal'])
        eta_e = np.matmul(Cov_Inv[label]['tumor'],mean[label]['tumor'])
        eta_c = eta_n + eta_e
        
        mixed = np.matmul(np.transpose(eta_n),np.matmul(covariance_matrix[label]['normal'],eta_n))\
               +np.matmul(np.transpose(eta_e),np.matmul(covariance_matrix[label]['tumor'],eta_e))\
               -np.matmul(np.transpose(eta_c),np.matmul(Vc,eta_c))\
        
        arg= -0.5*(DIMENSION*math.log(2*math.pi) + log_sum + mixed)
        
        
        overlap[label]= math.pow(2,A)*math.pow(2*math.pi,B)*math.exp(arg)*math.pow(det_inv_Vc,-0.25)
        log_overlap[label]=-A*math.log(2)-B*math.log(2*math.pi)-arg +0.25*math.log(det_inv_Vc)
    return overlap, log_overlap



def _complexity_map(mean, covariance_matrix ,DIMENSION, EXCLUSION_LIST):  
    entropy_dict = _H_multiv_gaussian(covariance_matrix,DIMENSION)
    
    delta_H = _delta_H(entropy_dict)
    overlap, log_overlap = _overlap(mean, covariance_matrix ,DIMENSION)
    
    x1,x2,x3,y,z= [],[],[],[],[]
    
    for label in log_overlap.keys():
        if (label in EXCLUSION_LIST):
            pass
        else:
            y.append(log_overlap[label])
            x1.append(entropy_dict[label]['normal'])
            x2.append(entropy_dict[label]['tumor'])
            x3.append(delta_H[label])
            z.append(label)
        
    
    fig, ax = plt.subplots()
    ax.scatter(x3, y, label=r'$\Delta S=S_{tumor}-S_{normal}$')
    ax.scatter(x1, y, label=r'$S_{normal}$')
    ax.scatter(x2, y, label=r'$S_{tumor}$')
    
    for i, txt in enumerate(z):
        ax.annotate(txt, (x1[i], y[i]))
        ax.annotate(txt, (x2[i], y[i]))
        ax.annotate(txt, (x3[i], y[i]))
    
    ident = [0.0, 100.0]
    fit= np.polyfit(x3,y,1)#Polynomial coefficients, highest power first
    
    fitted_y =[fit[0]*i+fit[1] for i in ident] 
    
    
    plt.plot(ident,fitted_y,'b--',label=r"$m \approx {0},n \approx {1}$".format(np.around(fit[0], 2),np.around(fit[1], 2)),linewidth=1)
    plt.xlabel(r"$\Delta S, S$ (nats)")
    plt.ylim(0, 100)
    plt.ylabel(r"$-\ln \,I$")
    plt.title("Overlap $vs$ Entropy")
    pylab.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig('overlapentropy_fig.pdf')   
    
    return fit

# The main program
if __name__ == "__main__":
    data, missmatch_pc_samples=read_data(srcpc,srcsampl)
    core_matrix = _transpose_data(data)

    EXCLUSION_LIST =['READ','ESCA','CESC','GBM','PAAD','BLCA']# BECAUSE FEW 'normal' SAMPLES, OTHERWISE EMPTY. 

    _PC_pair_plot(core_matrix, 1, 2,'LUSC') 

    low_DIMENSION_core_matrix=_low_dimension_core_matrix(core_matrix, LOW_DIMENSION=20)
    low_DIMENSION_covariance_matrix = _covariance_matrix_(low_DIMENSION_core_matrix)
    low_DIMENSION_mean = _mean_dict(low_DIMENSION_core_matrix)
    low_DIMENSION_std = _std_dict(low_DIMENSION_core_matrix)

    linear_fit =_complexity_map(low_DIMENSION_mean,  low_DIMENSION_covariance_matrix ,20, EXCLUSION_LIST)

