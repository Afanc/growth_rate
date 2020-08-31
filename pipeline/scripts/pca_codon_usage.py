#!/home/dmollet/.conda/envs/python_env/bin/python3

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="trna analysis"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

import glob
import os
import sys
import csv
import operator
import numpy as np
from sklearn.decomposition import FactorAnalysis, PCA
from sklearn.preprocessing import StandardScaler
from Bio import SeqIO
from Bio.SeqUtils import GC_skew
from fastdtw import fastdtw
from dtw import dtw
from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans
from skbio.stats.ordination import pcoa
import scipy.cluster.hierarchy as spc
import matplotlib.pyplot as plt
from dict_of_codons import * 

#I know I'll divide by zero thanks
#np.seterr(divide='ignore', invalid='ignore')

trans = str.maketrans('ATGC', 'TACG')

# input files 
file_list = glob.glob('fasta/*.fasta')
file_list = sorted(file_list)
print(file_list)

# get data
data = []
names = []
for file_path in file_list:
    for record in SeqIO.parse(file_path, "fasta"):
        names.append(record.id)
        data.append(record.seq)

number_strains = len(data)

codon_usages = []
dicodon_usages = []

for num,i in enumerate(data) :

    temp = CodonsDict.copy()
    ditemp = diCodonsDict.copy()

    str_seq = str(i)
    
    f = "stats_all"+file_list[num][5:-5]+"stats_all"
    
    print(f)
    gene_positions = []

    with open(f) as s:
        stat_file = csv.reader(s, delimiter=' ')
        next(stat_file, None)         #header

        for j in stat_file :
            if(j[1] != '0') :
                gene_positions.append([int(j[2]), int(j[3]), int(j[4])])

    codons = []
    dicodons = []
    summ = 0
    #print(gene_positions)

    for k in gene_positions :

        summ += int(k[1]) - int(k[0]) + 1

        if k[2] == 1 :
            codons.extend([str_seq[j:j+3] for j in range(k[0]-1, k[1]-2, 3)]) #range(0,len(str_seq) - len(str_seq)%3, 3)]
            dicodons.extend([str_seq[j:j+6] for j in range(k[0]-1, k[1]-5, 3)]) #range(0,len(str_seq) - len(str_seq)%3, 3)]

        elif k[2] == 0:
            codons.extend([str_seq[j-3:j][::-1].translate(trans) for j in range(k[1], k[0]+1, -3)]) #range(0,len(str_seq) - len(str_seq)%3, 3)]
            dicodons.extend([str_seq[j-6:j][::-1].translate(trans) for j in range(k[1], k[0]+4, -3)]) #range(0,len(str_seq) - len(str_seq)%3, 3)]
            
    #print(codons)

    for codon in codons :
        c = codon.upper()
        try: 
            temp[c]+=1
        except:
            pass
            #print("codon not in dict : ", c)

    codon_usages.append(temp)

    for dicodon in dicodons :
        c = dicodon.upper()
        try: 
            ditemp[c]+=1
        except:
            pass

    dicodon_usages.append(ditemp)

    print("codon usage sum", sum(codon_usages[num].values()))
    print("len of genome", summ/3)
   

conc_data = np.zeros(shape=(len(codon_usages), 64))
diconc_data = np.zeros(shape=(len(dicodon_usages), 4096))

for i,dist in enumerate(codon_usages) :
    conc_data[i,] = np.fromiter(dist.values(), dtype=float)

for i,dist in enumerate(dicodon_usages) :
    diconc_data[i,] = np.fromiter(dist.values(), dtype=float)

#PCA
pca = PCA(svd_solver = 'full')
scaler = StandardScaler()

normed_data = scaler.fit_transform(conc_data)

res_pca = pca.fit_transform(normed_data)

#---------------------------------K-Means final
kmeans = KMeans(init='k-means++', n_clusters=3, n_init=1000)
res_kmeans = kmeans.fit(normed_data)

#---------------------------------Plottig

fig = plt.figure()
y_pos = np.arange(len(pca.explained_variance_ratio_))
plt.bar(y_pos, pca.explained_variance_ratio_)
fig.suptitle('PC Explained Variance [codon usage]', fontsize=16)
fig.savefig('results/explained_var_codon_usage.png')

fig, ax = plt.subplots(3, figsize = (8,16))
ax[0].scatter(res_pca[:,0], res_pca[:,1], c=res_kmeans.labels_)
ax[1].scatter(res_pca[:,1], res_pca[:,2], c=res_kmeans.labels_)
ax[2].scatter(res_pca[:,2], res_pca[:,3], c=res_kmeans.labels_)
ax[0].set(xlabel='PC 1, '+str(round(100*pca.explained_variance_ratio_[0],1))+'% var', ylabel='PC 2, '+str(round(100*pca.explained_variance_ratio_[1],1))+'% var')
ax[1].set(xlabel='PC 2, '+str(round(100*pca.explained_variance_ratio_[1],1))+'% var', ylabel='PC 3, '+str(round(100*pca.explained_variance_ratio_[2],1))+'% var')
ax[2].set(xlabel='PC 3, '+str(round(100*pca.explained_variance_ratio_[2],1))+'% var', ylabel='PC 4, '+str(round(100*pca.explained_variance_ratio_[3],1))+'% var')

for i, txt in enumerate(names):
    ax[0].annotate(txt, (res_pca[i,0], res_pca[i,1]))
for i, txt in enumerate(names):
    ax[1].annotate(txt, (res_pca[i,1], res_pca[i,2]))
for i, txt in enumerate(names):
    ax[2].annotate(txt, (res_pca[i,2], res_pca[i,3]))

fig.suptitle('PCA - [codon usage]', fontsize=16)

fig.savefig('results/pca_codon_usage.png')

#PCA
pca = PCA(svd_solver = 'full')
scaler = StandardScaler()

normed_data = scaler.fit_transform(diconc_data)

res_pca = pca.fit_transform(normed_data)

#---------------------------------K-Means final
kmeans = KMeans(init='k-means++', n_clusters=3, n_init=1000)
res_kmeans = kmeans.fit(normed_data)

#---------------------------------Plottig

fig = plt.figure()
y_pos = np.arange(len(pca.explained_variance_ratio_))
plt.bar(y_pos, pca.explained_variance_ratio_)
fig.suptitle('PC Explained Variance [dicodon usage]', fontsize=16)
fig.savefig('results/explained_var_dicodon_usage.png')

fig, ax = plt.subplots(3, figsize = (8,16))
ax[0].scatter(res_pca[:,0], res_pca[:,1], c=res_kmeans.labels_)
ax[1].scatter(res_pca[:,1], res_pca[:,2], c=res_kmeans.labels_)
ax[2].scatter(res_pca[:,2], res_pca[:,3], c=res_kmeans.labels_)
ax[0].set(xlabel='PC 1, '+str(round(100*pca.explained_variance_ratio_[0],1))+'% var', ylabel='PC 2, '+str(round(100*pca.explained_variance_ratio_[1],1))+'% var')
ax[1].set(xlabel='PC 2, '+str(round(100*pca.explained_variance_ratio_[1],1))+'% var', ylabel='PC 3, '+str(round(100*pca.explained_variance_ratio_[2],1))+'% var')
ax[2].set(xlabel='PC 3, '+str(round(100*pca.explained_variance_ratio_[2],1))+'% var', ylabel='PC 4, '+str(round(100*pca.explained_variance_ratio_[3],1))+'% var')

for i, txt in enumerate(names):
    ax[0].annotate(txt, (res_pca[i,0], res_pca[i,1]))
for i, txt in enumerate(names):
    ax[1].annotate(txt, (res_pca[i,1], res_pca[i,2]))
for i, txt in enumerate(names):
    ax[2].annotate(txt, (res_pca[i,2], res_pca[i,3]))

fig.suptitle('PCA - [dicodon usage]', fontsize=16)

fig.savefig('results/pca_dicodon_usage.png')
