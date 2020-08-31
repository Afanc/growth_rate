#!/home/dmollet/.conda/envs/python_env/bin/python3 -u

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="GC_skew"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

import glob
import os
import sys
import operator
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC_skew
from fastdtw import fastdtw
from dtw import dtw
from scipy.spatial.distance import euclidean
from sklearn.cluster import KMeans
from skbio.stats.ordination import pcoa
import scipy.cluster.hierarchy as spc
import matplotlib.pyplot as plt

#I know I'll divide by zero thanks
#np.seterr(divide='ignore', invalid='ignore')

# input files 
file_list = glob.glob('fasta/*.fasta')
file_list = sorted(file_list)

# get data
data = []
names = []
for file_path in file_list:
    for record in SeqIO.parse(file_path, "fasta"):
        names.append(record.id)
        data.append(record.seq)

print(names)
number_strains = len(data)

GC_skews = []
for i in data :
    GC_skews.append(GC_skew(i, window=200))


#comparing specific strains
if True:

    fig = plt.figure(figsize=(10,10))

    pos_jvci1 = names.index("JCVI1")
    pos_jvci2 = names.index("JCVI2")
    pos_jvci3 = names.index("JCVI3")
    pos_gm12wt = names.index("GM12WT")
    pos_gm12d86 = names.index("GM12delta68")
    pos_mg37 = names.index("MG37")
    pos_Mf5583 = names.index("Mf5583")

    compared = [pos_jvci1, pos_jvci2, pos_jvci3, pos_gm12wt, pos_gm12d86, pos_mg37, pos_Mf5583]
    max_length = np.max([len(GC_skews[i]) for i in compared])

    smallest_seq = np.argmin([len(str(data[i])) for i in compared])

    for i in compared :
        cum_sum = np.cumsum(GC_skews[i])

        #two next lines : interpolation
        #interpolated = interp.interp1d(np.arange(cum_sum.size),cum_sum)
        #cum_sum = interpolated(np.linspace(0,cum_sum.size-1,max_length))

        #next lines, adaptive window GC_skew - NOT GOOD
        #adapted_w = 200*len(str(data[i]))//len(str(data[smallest_seq]))
        #GC_skews_adaptive = GC_skew(data[i], window=adapted_w)
        #cum_sum = np.cumsum(GC_skews_adaptive)
        #print(len(cum_sum)) #ACH not the same length, be more precise

        plt.plot(cum_sum, label=names[i])

    leg = plt.legend(prop={'size':16})

    for line in leg.get_lines():
        line.set_linewidth(4.0)

    plt.xlabel('bp', fontsize=16)
    plt.ylabel('CGCSkew', fontsize=16)

    #fig.suptitle('Cumulative GC Skew', fontsize=16)
    fig.savefig('results/GCC.png')

    #plt.show()


all_distances = np.zeros(shape=(number_strains, number_strains))
l = np.empty(number_strains)

for key_1, strain_1 in enumerate(GC_skews) :
    l[key_1] = len(GC_skews[key_1])
    for key_2 in range(key_1+1, number_strains) :
        dist, path = fastdtw(GC_skews[key_1], GC_skews[key_2])

        #print("dist", dist)
        all_distances[key_1, key_2] = dist#/(len(GC_skews[key_1]) + len(GC_skews[key_2]))
        all_distances[key_2, key_1] = dist#/(len(GC_skews[key_1]) + len(GC_skews[key_2]))

print(all_distances.shape)


#PCoA
pcoa_res = pcoa(all_distances, method='eigh')
print(pcoa_res)

pcoa_samples = pcoa_res.samples.to_numpy()

fig, ax = plt.subplots(3, figsize = (12,24))

ax[0].scatter(pcoa_samples[:,0], pcoa_samples[:,1])
ax[1].scatter(pcoa_samples[:,1], pcoa_samples[:,2])
ax[2].scatter(pcoa_samples[:,2], pcoa_samples[:,3])

ax[0].set(xlabel='PC 1, '+str(round(100*pcoa_res.proportion_explained[0],1))+'% var', ylabel='PC 2, '+str(round(100*pcoa_res.proportion_explained[1],1))+'% var')
ax[1].set(xlabel='PC 2, '+str(round(100*pcoa_res.proportion_explained[1],1))+'% var', ylabel='PC 3, '+str(round(100*pcoa_res.proportion_explained[2],1))+'% var')
ax[2].set(xlabel='PC 3, '+str(round(100*pcoa_res.proportion_explained[2],1))+'% var', ylabel='PC 4, '+str(round(100*pcoa_res.proportion_explained[3],1))+'% var')

for i, txt in enumerate(names):
    ax[0].annotate(txt, (pcoa_samples[i,0], pcoa_samples[i,1]))
for i, txt in enumerate(names):
    ax[1].annotate(txt, (pcoa_samples[i,1], pcoa_samples[i,2]))
for i, txt in enumerate(names):
    ax[2].annotate(txt, (pcoa_samples[i,2], pcoa_samples[i,3]))

#fig.suptitle('PCoA on DTW of GCSkews', fontsize=16)

fig.savefig('results/pcoa_gcskew.png')

# cluster
pdist = spc.distance.pdist(all_distances) #pairwise distances

print("pdist", pdist)

linkage = spc.linkage(pdist, method='average') 
print("linkage", linkage)

print("\n")
s = np.sum(all_distances,0)
print(s, l)

print(np.corrcoef(s, l))

#exit()
fig = plt.figure(figsize=(20,20))

spc.dendrogram(linkage, labels=names)

fig.suptitle('DTW of GCSkew')
fig.savefig('results/DTW_GCSkew.png')

fig, ax = plt.subplots(figsize = (12,12))

ax.scatter(pcoa_samples[:,0], pcoa_samples[:,1], s=400)
ax.set_xlabel('PC 1, '+str(round(100*pcoa_res.proportion_explained[0],1))+'% var', fontsize=18)
ax.set_ylabel('PC 2, '+str(round(100*pcoa_res.proportion_explained[1],1))+'% var', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.xaxis.set_label_coords(0.5, -0.07)

for i, txt in enumerate(names):
    ax.annotate(txt, (pcoa_samples[i,0]+1.5, pcoa_samples[i,1]), fontsize=24)

fig.savefig('results/pcoa_gcskew.png')


exit()

#cor = df.corr()
#cor = cor.dropna(axis=0, how='all').dropna(axis=1, how='all')

#d, cost_matrix, acc_cost_matrix, path = dtw(a, b, dist=euclidean)
#print("dtw", d)
