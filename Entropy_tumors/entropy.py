"""
This script file computes the entropy and overlap of PCA. 
The procedure and main results are described in paper 
arXiv:2005.02271
For more details on this file see author(s):
FEQ, DALV
"""

import numpy as np
import math, pylab
import xlrd as xl
import matplotlib.pyplot as plt
import sys,os
import random

# General variables ------------------------------------------------------

Npc = 20
tissues_id = ['BRCA','COAD','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC']
N_random = 1000

srcpc = "../databases_generated/TCGA_pca/"
srcsampl = "../databases_external/TCGA/"
sample_path = '../databases_external/TCGA/'
PC_path = '../databases_generated/TCGA_pca/'


# Functions ------------------------------------------------

def read_PC(sample_path,PC_path,tissue_id,Npc):
    pc_normal = []
    pc_tumor = []
    sample_wb = xl.open_workbook(sample_path+"sample"+tissue_id+".xls")
    sample_ws = sample_wb.sheet_by_index(0)
    pc_wb = xl.open_workbook(PC_path+"pc"+tissue_id+".xls")
    pc_ws = pc_wb.sheet_by_index(0)
    print('Total of',tissue_id+' cases:', sample_ws.nrows)
    for i in range(sample_ws.nrows-1):
        if(sample_ws.cell_value(i+1,3) == 'Solid Tissue Normal'):
            pc_normal.append([])
            for j in range(Npc):    
                pc_normal[-1].append(pc_ws.cell_value(i,j))
        else:
            pc_tumor.append([])
            for j in range(Npc):    
                pc_tumor[-1].append(pc_ws.cell_value(i,j))

    pc_normal = np.array(pc_normal, dtype=np.float64)
    pc_tumor  = np.array(pc_tumor , dtype=np.float64)

    return pc_normal, pc_tumor


def compute_entropy(pc_dataset):

    coreM=pc_dataset.T
    cov_mat=np.cov(coreM,ddof=1)
    mean=np.mean(coreM,axis=1)
    stdev=np.std(coreM,axis=1)
    u,ss,vh =np.linalg.svd(cov_mat)
    s=np.diag(ss)
    cov_inv=np.matmul(np.matmul(vh.T,np.linalg.inv(s)),u.T)
    sum_log=sum(math.log(s[i][i]) for i in range(Npc))
    S_dataset=0.5*(sum_log + Npc*(1 + math.log(2*math.pi)))

    return S_dataset,mean,cov_mat,cov_inv


def compute_overlap(Npc,mean_n,mean_t,cov_mat_n,cov_mat_t,cov_inv_n,cov_inv_t):

    A, B = Npc/2 , Npc/4
    inv_Vc = cov_inv_n+cov_inv_t
    Vc = np.linalg.inv(inv_Vc)
        
    u_Vc, s_Vc, vh_Vc = np.linalg.svd(Vc)
    u_inv_Vc, s_inv_Vc, vh_inv_Vc = np.linalg.svd(inv_Vc)
    det_Vc, det_inv_Vc = np.prod(s_Vc), np.prod(s_inv_Vc)
        
    u_inv_Vn, s_inv_Vn, vh_inv_Vn = np.linalg.svd(cov_inv_n)
    u_inv_Ve, s_inv_Ve, vh_inv_Ve = np.linalg.svd(cov_inv_t)        
        
    log_sum=sum(-math.log(s_inv_Vn[i])-math.log(s_inv_Ve[i])+math.log(s_inv_Vc[i]) for i in range(Npc))        
        
    eta_n = np.matmul(cov_inv_n,mean_n)
    eta_e = np.matmul(cov_inv_t,mean_t)
    eta_c = eta_n + eta_e
        
    mixed = np.matmul(np.transpose(eta_n),np.matmul(cov_mat_n,eta_n))\
           +np.matmul(np.transpose(eta_e),np.matmul(cov_mat_t,eta_e))\
           -np.matmul(np.transpose(eta_c),np.matmul(Vc,eta_c))\
        
    arg= -0.25*(Npc*math.log(2*math.pi) + log_sum + mixed)        
        
    overlap = math.pow(2,A)*math.pow(2*math.pi,B)*math.exp(arg)*math.pow(det_inv_Vc,-0.25)
    log_overlap = -A*math.log(2)-B*math.log(2*math.pi)-arg +0.25*math.log(det_inv_Vc)

    return overlap, log_overlap


# The main program -----------------------------------------------------------------------------------------
if __name__ == "__main__":
    random.seed(a=1234567890)
    #random.seed()
    N_normal,N_tumor=[],[]
    S_normal,S_tumor,S_delta=[],[],[]
    I,lnI=[],[]
    S_normal_samp_mean,S_normal_samp_std,S_tumor_samp_mean,S_tumor_samp_std,S_delta_samp=[],[],[],[],[]
    S_tumor_samp_min_mean,S_tumor_samp_min_std,S_delta_samp_min=[],[],[]
    lnI_samp_mean,lnI_samp_std,lnI_samp_min_mean,lnI_samp_min_std=[],[],[],[]
    Nn_min=23
    print('\nLoading TCGA Principal Components\n')
    for t_id in tissues_id: 
        pc_normal,pc_tumor=read_PC(sample_path,PC_path,t_id,Npc)
        N_normal.append(len(pc_normal))
        N_tumor.append(len(pc_tumor))
        print('normal/tumor:',len(pc_normal),'/',len(pc_tumor))
        print('PCs:',str(len(pc_normal[0]))+'x'+str(len(pc_tumor[0])))

        Sn_i,St_i,lnI_i=[],[],[]
        for i in range(N_random):
            pc_normal_samp=pc_normal[random.sample(range(len(pc_normal)),Nn_min)]
            pc_tumor_samp=pc_tumor[random.sample(range(len(pc_tumor)),len(pc_normal))]
            S_n,mean_n,cov_mat_n,cov_inv_n=compute_entropy(pc_normal_samp)
            S_t,mean_t,cov_mat_t,cov_inv_t=compute_entropy(pc_tumor_samp)
            Sn_i.append(S_n)
            St_i.append(S_t)
            I_id,logI_id = compute_overlap(Npc,mean_n,mean_t,cov_mat_n,cov_mat_t,cov_inv_n,cov_inv_t)
            lnI_i.append(logI_id)
        S_normal_samp_mean.append(np.mean(Sn_i))
        S_normal_samp_std.append(np.std(Sn_i))
        S_tumor_samp_min_mean.append(np.mean(St_i))
        S_tumor_samp_min_std.append(np.std(St_i))
        lnI_samp_min_mean.append(np.mean(lnI_i))
        lnI_samp_min_std.append(np.std(lnI_i))


        S_n,mean_n,cov_mat_n,cov_inv_n=compute_entropy(pc_normal)
        S_normal.append(S_n)
        S_t,mean_t,cov_mat_t,cov_inv_t=compute_entropy(pc_tumor)
        S_tumor.append(S_t)
        I_id,logI_id = compute_overlap(Npc,mean_n,mean_t,cov_mat_n,cov_mat_t,cov_inv_n,cov_inv_t)
        I.append(I_id)
        lnI.append(logI_id)

        S_i,lnI_i=[],[]
        for i in range(N_random):
            pc_tumor_samp=pc_tumor[random.sample(range(len(pc_tumor)),len(pc_normal))]
            S_t,mean_t,cov_mat_t,cov_inv_t=compute_entropy(pc_tumor_samp)
            S_i.append(S_t)
            I_id,logI_id = compute_overlap(Npc,mean_n,mean_t,cov_mat_n,cov_mat_t,cov_inv_n,cov_inv_t)
            lnI_i.append(logI_id)
        S_tumor_samp_mean.append(np.mean(S_i))
        S_tumor_samp_std.append(np.std(S_i))
        lnI_samp_mean.append(np.mean(lnI_i))
        lnI_samp_std.append(np.std(lnI_i))

        S_delta.append(S_tumor[-1]-S_normal[-1])
        S_delta_samp.append(S_tumor_samp_mean[-1]-S_normal[-1])
        S_delta_samp_min.append(S_tumor_samp_min_mean[-1]-S_normal_samp_mean[-1])

    np.savetxt('entropies_dat.dat', np.transpose([tissues_id,N_normal,N_tumor,S_normal,S_tumor,S_delta,
                                    I,lnI]),
    header='tissue\tNnormal\tNtumor\tSn\t\t\tSt\t\t\tdeltaS\t\t\tI\t\t\t-ln(I)',
    comments='', delimiter=" \t", fmt="%s")
    np.savetxt('entropies_samp-'+str(Nn_min)+'_dat.dat', np.transpose([tissues_id,N_normal,N_tumor,
    S_normal_samp_mean,S_normal_samp_std,S_tumor_samp_min_mean,S_tumor_samp_min_std,S_delta_samp_min,
    lnI_samp_min_mean,lnI_samp_min_std]),
    header='tissue\tNnormal\tNtumor\t<Sn>_'+str(Nn_min)+'\t\t\tstd(Sn)_'+str(Nn_min)+'\t\t<St>_'+
    str(Nn_min)+'\t\t\tstd(St)_'+str(Nn_min)+'\t\t<St>-<Sn>\t\t-<lnI>_'+str(Nn_min)+'\t\tstd(lnI)_'+
    str(Nn_min),comments='',delimiter=" \t",fmt="%s")
    np.savetxt('entropies_samp-Nn_dat.dat', np.transpose([tissues_id,N_normal,N_tumor,S_normal,
    S_tumor_samp_mean,S_tumor_samp_std,S_delta_samp,lnI_samp_mean,lnI_samp_std]),
    header='tissue\tNnormal\tNtumor\tSn\t\t\t<St>_Nn\t\t\tstd(St)_Nn\t\t<St>-Sn\t\t\t-<lnI>_Nn'+
    '\t\tstd(lnI)_Nn',comments='',delimiter=" \t",fmt="%s")


    print('\nPloting figures')
# plotting PC1 vs. PC2 ----------------------------
    t_id='LUSC'
    pc_normal,pc_tumor=read_PC(sample_path,PC_path,t_id,2)
    plt.scatter(pc_normal[:,0],pc_normal[:,1],label='normal',c='b',s=15)
    plt.scatter(pc_tumor[:,0],pc_tumor[:,1],label='tumor',c='r',s=15)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title(t_id)
    plt.legend()
    plt.tight_layout()
    plt.grid(True)
    plt.savefig('ge'+t_id+'_pc1-2_fig.pdf')
    plt.clf()


# plotting -lnI vs. deltaS ----------------------------
    fit=np.polyfit(S_normal,lnI_samp_mean,1)
    f=np.poly1d(fit)
    plt.errorbar(S_normal,lnI_samp_mean,yerr=lnI_samp_std,fmt='o',
    capsize=5,label='normal',c='blue')
    plt.errorbar(S_tumor_samp_mean,lnI_samp_mean,xerr=S_tumor_samp_std,yerr=lnI_samp_std,fmt='o',
    capsize=5,label='tumor',c='red')
    plt.plot([min(S_normal),max(S_normal)],f([min(S_normal),max(S_normal)]),c='blue',linestyle='dashed')
    plt.axvline(x=np.mean(S_tumor_samp_mean),color='red',linestyle='dashed')
    plt.xlabel('Sn,<St>')
    plt.ylabel('-lnI')
    plt.legend()
    plt.tight_layout()
    plt.legend(loc='upper center')
    plt.savefig('entropies_map_fig.pdf')
    plt.clf()

    fit=np.polyfit(S_delta_samp,lnI_samp_mean,1)
    f=np.poly1d(fit)
    plt.errorbar(S_delta_samp,lnI_samp_mean,xerr=S_tumor_samp_std,yerr=lnI_samp_std,fmt='o',
    capsize=5,label='Nn random samples',c='black')
    plt.plot([min(S_delta_samp),max(S_delta_samp)],f([min(S_delta_samp),max(S_delta_samp)]),
    c='red',linestyle='dashed')
    for i in range(len(S_delta_samp)):
        plt.annotate(tissues_id[i], (S_delta_samp[i], lnI_samp_mean[i]), ha='center')
    plt.xlabel('<St>-Sn')
    plt.ylabel('-lnI')
    plt.legend()
    plt.tight_layout()
    plt.savefig('entropy_samp_map_fig.pdf')
    plt.clf()

print("Done!\n")
print("Check out the outputs!\n")
