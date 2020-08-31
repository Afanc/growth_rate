#!/home/dmollet/.conda/envs/python_env/bin/python3

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="correlations"
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
import scipy.cluster.hierarchy as spc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#I know I'll divide by zero thanks
np.seterr(divide='ignore', invalid='ignore')

# input files 
file_list = glob.glob('stats/*.stats')
file_list = sorted(file_list)

# get stuff from input
header_dic = {'rank':1, 'beg':2, 'end':3, 'dir':4, 'occ':5}
header_col = tuple(header_dic[i] for i in sys.argv[1:])
params = sys.argv[1:]
if not header_col :
    header_col = (2,4,5) #default is rank, dir, occ
    params = 'beg,dir,occ'

# get data
data = []
names = []
for file_path in file_list:
    names.append(os.path.splitext(file_path)[0][6:])
    data.append(np.genfromtxt(file_path, delimiter=' ', skip_header=1, usecols=header_col, case_sensitive=True, dtype=None))

# all the gene_names in one nice vector
gene_names = np.array(np.genfromtxt(file_list[0], delimiter=' ', skip_header=1, usecols=(0), case_sensitive=True, dtype='str'), dtype='object')

# replication rates
repl_rates = np.array(np.genfromtxt('repl_rates.csv', delimiter=',', skip_header=0, case_sensitive=True, dtype='str'), dtype='object')

# foolproofing
if (repl_rates.shape[0] != len(names)):
    print("oooooh shit ! \n# replication rates != # strains !\n")
    exit()

# sort replication rates too
repl_rates = repl_rates[np.argsort(repl_rates[:, 0])]

# and get all the features in one hell of a vector
num_genes = gene_names.shape[0]

features = np.empty((len(header_col)*num_genes), dtype="object")

i = 0
for key, value in sorted(header_dic.items(), key=operator.itemgetter(1)):
    if value in header_col :
        features[num_genes*(i):num_genes*(i+1)] = gene_names+"_"+key
        i+=1

# ------ all the data is a single vector

data = np.asarray(data, dtype=np.float32)
conc_data = np.concatenate(np.swapaxes(data,0,2), axis=0).transpose()

df = pd.DataFrame(conc_data, columns = features)

# ---- simple correlation

#normalize
scaler = StandardScaler()
repl_rates = scaler.fit_transform(repl_rates[:,1].reshape(-1,1))
conc_data = scaler.fit_transform(conc_data) #doesn't normalize sample-wise

#normalize again ? somehow standardscaler is not sufficient
#repl_m = repl_rates[:,1][:,None].T.astype(np.float) #if not using the standardscaler
repl_m = repl_rates.T.astype(np.float)
conc_m = (conc_data.T - conc_data.mean(0)[:,None])
#conc_m = (conc_data - conc_data.mean(1)[:,None]).T #before, weird results ?

# Sum of squares across rows
ss_repl = (repl_m**2).sum(1)
ss_conc = (conc_m**2).sum(1)

# Finally get corr coeff
features_cor = np.dot(repl_m, conc_m.T)/np.sqrt(np.dot(ss_repl[:,None],ss_conc[None]))

# Sort
sorted_features_idx = np.argsort(features_cor)

sorted_features = features[sorted_features_idx][0]
sorted_cor = features_cor[0][sorted_features_idx][0]

#remove NAs
sorted_features = sorted_features[~np.isnan(sorted_cor)]
sorted_cor = sorted_cor[~np.isnan(sorted_cor)]

# plot
fig, ax = plt.subplots(1, figsize = (40,40))

y_pos = np.arange(sorted_features.shape[0])
plt.barh(y_pos, sorted_cor, tick_label=sorted_features)
ax.tick_params(axis='both', which='major', labelsize=8)

fig.suptitle('Correlations - features : '+str(params), fontsize=36)
fig.savefig("results/features_correlations.png")

print(sorted_features[:40])
print(sorted_cor[:40])
print(sorted_features[:-40:-1])
print(sorted_cor[:-40:-1])

exit()

# ---- multiple correlations and hierarchical clustering
# correlation first
cor = df.corr()

cor = cor.dropna(axis=0, how='all').dropna(axis=1, how='all')

# cluster
pdist = spc.distance.pdist(cor) #pairwise distances
linkage = spc.linkage(pdist, method='average') #clustering. Using average distance cluster vs cluster (non-weighted by number of elements)
idx = spc.fcluster(linkage, 0.5*pdist.max(), 'distance')  #how to separate clusters ? #distance : 0.5 * pdist.max() #maxclust : 3

# get indexes out of clusters
columns = [df.columns.tolist()[i] for i in list((np.argsort(idx)))] #just the columns' names in the right/clustered order

# and reorder
df = df.reindex(columns, axis=1) #reorder df by column order. Not a square matrix

# before doing correlation again...
cor = df.corr()
cor = cor.dropna(axis=0, how='all').dropna(axis=1, how='all')

# plot this shit
fig = plt.figure(figsize=(120,100))                                                                                                                                                                         
sns.heatmap(cor, annot=False, cmap=plt.cm.Reds)

fig.suptitle('Cross-correlations - features : '+str(params), fontsize=36)
fig.savefig('results/correlated_features.png')                                                                                                                                                                      
