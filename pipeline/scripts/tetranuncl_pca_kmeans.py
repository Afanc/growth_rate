#!/home/dmollet/.conda/envs/python_env/bin/python3 -u

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="pca/km"
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
from sklearn.decomposition import FactorAnalysis, PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from pandas import *
import matplotlib.pyplot as plt

#get data
file_list = glob.glob('tetranucl/*.tetranucl')
file_list = sorted(file_list)

data = []
names = []
for file_path in file_list:
    names.append(os.path.splitext(file_path)[0][10:])
    data.append(np.genfromtxt(file_path, delimiter=',', skip_header=1, usecols=(1), case_sensitive=True, dtype=None))

data = np.asarray(data, dtype=np.float32)

#conc_data = np.concatenate(np.swapaxes(data,0,2), axis=0).transpose()

# all the features in one nice vector
tetra_names = np.array(np.genfromtxt(file_list[0], delimiter=',', skip_header=1, usecols=(0), case_sensitive=True, dtype='str'), dtype='object')

#---------------------------------PCA
pca = PCA(svd_solver = 'full')
scaler = StandardScaler()

normed_data = scaler.fit_transform(data)

res_pca = pca.fit_transform(normed_data)

# representative features

#first_pc = np.transpose(np.vstack([features, pca.components_[0,...]]))
#first_pc = first_pc[first_pc[...,1].argsort(axis=0)]
#
#second_pc = np.transpose(np.vstack([features, pca.components_[1,...]]))
#second_pc = second_pc[second_pc[...,1].argsort(axis=0)]

#np.savetxt('results/first_pc.csv', X=first_pc, fmt=['%s', '%.4e'], delimiter=',', header="feature, PC1 : "+str(params))
#np.savetxt('results/second_pc.csv', X=second_pc, fmt=['%s', '%.4e'], delimiter=',', header="feature, PC2 : "+str(params))

#---------------------------------K-Means
if False :
    fig = plt.figure()
    kmeans_list = []
    clusters = [2,3,4,5,6,7]
    cluster_names = ['  '+str(i)+'  ' for i in clusters] #STUPID

    for i in clusters :
        kmeans = KMeans(init='k-means++', n_clusters=i, n_init=1000)
        kmeans_list.append((kmeans.fit(normed_data).labels_).tolist())
        #colors instead of labels (for plotting)
        kmeans_list[-1] = ['r' if i == 0 else 'b' if i == 1 else 'g' if i == 2 else 'y' if i == 3 else 'w' if i == 4 else 'c' if i == 5 else 'k' for i in kmeans_list[-1]]

    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111, frameon=False, xticks=[], yticks=[])
    the_table = plt.table(cellText=None, rowLabels=cluster_names, colLabels=names, cellColours=kmeans_list, loc='center')

    # ALL THE STUPID FORMATTING
    the_table.scale(1.2,3.5)
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(9)

    for cell in the_table._cells:
        if cell[0] ==0:
            the_table._cells[cell].get_text().set_rotation(45)
    # ACH

    fig.suptitle('K-means - features : '+str(params), fontsize=16)
    fig.savefig('results/kmeans.png')

#---------------------------------K-Means final
kmeans = KMeans(init='k-means++', n_clusters=4, n_init=1000)
res_kmeans = kmeans.fit(normed_data)

#---------------------------------Plottig

fig = plt.figure()
y_pos = np.arange(len(pca.explained_variance_ratio_))
plt.bar(y_pos, pca.explained_variance_ratio_)
fig.suptitle('PC Explained Variance - Tetranucleotides', fontsize=16)
fig.savefig('results/tetranucl_explained_var.png')

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

fig.suptitle('PCA - Tetranucleotides', fontsize=16)

fig.savefig('results/tetranucl_pca_1_2_3_4.png')

fig, ax = plt.subplots(figsize = (12,12))
ax.scatter(res_pca[:,0], res_pca[:,1], c=res_kmeans.labels_, s=400)
ax.set_xlabel('PC 1, '+str(round(100*pca.explained_variance_ratio_[0],1))+'% var', fontsize=18)
ax.set_ylabel('PC 2, '+str(round(100*pca.explained_variance_ratio_[1],1))+'% var', fontsize=18)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.xaxis.set_label_coords(0.5, -0.07)


for i, txt in enumerate(names):
    ax.annotate(txt, (res_pca[i,0]+1.5, res_pca[i,1]), fontsize=24)

#fig.suptitle('PCA - Tetranucleotides', fontsize=32)

fig.savefig('results/tetranucl_pca_1.png')


