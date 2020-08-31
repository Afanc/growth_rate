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
import re
import operator
import numpy as np

# input files 
file_list = glob.glob('full_genome/*.gen')
file_list = sorted(file_list)

# get data
data = []
for file_path in file_list:
    data.append(np.array(np.genfromtxt(file_path, delimiter=' ', skip_header=1, usecols=(0), case_sensitive=True, dtype='str'), dtype='object'))

data = [[re.sub('_\d+$', '', i) for i in j] for j in data] #remove indication of multiple occurrence
joined_genes = list(set().union(*data))

#print(joined_genes)

out_file = open("all_genes.txt","w")
out_file.writelines([i+"\n" for i in joined_genes])
