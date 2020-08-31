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
from Bio import SeqIO
from get_codon_usage_wide import *

#I know I'll divide by zero thanks
#np.seterr(divide='ignore', invalid='ignore')

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

print(file_list)

all_images = get_codon_usage(data, file_list)

#get replication rates (sort em !)
a = np.loadtxt('repl_rates.csv', dtype=np.str, delimiter=',')
a = sorted(a, key=lambda x: x[0])
repl_rates = [int(i[1]) for i in a]

np.savez('test', images=all_images, rates=repl_rates)

exit()
#to get visuals
for num,i in enumerate(all_images) :
    
    final_image = Image.fromarray(i).convert('L')         
    final_image.save('images/wide_'+file_list[num][6:-5]+'png')

    print("done !")

    sns.set_style('darkgrid')
    fig, axes = plt.subplots(nrows=1, ncols=len(all_images))
    fig.tight_layout()

    for i in range(0,len(all_images)):
        ax = fig.add_subplot(1,len(all_images),i+1)
        plt.imshow(all_images[i])#, cmap="Greys")
        plt.axis('off')
        axes[i].axis('off')
        ax.set_title(names[i], fontsize=8)

    fig.savefig('all_strains_wide.png', bbox_inches='tight')
