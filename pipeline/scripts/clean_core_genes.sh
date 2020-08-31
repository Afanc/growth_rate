#!/bin/bash
#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="clean"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

cat joined_core | sed 's/gene=//' > temp;
cat temp | grep -Fv "`cat non_core_genes.txt | cut -d" " -f3 | cut -d"=" -f2`" > joined_core ; 
sed -i -e 's/^/gene=/' joined_core;
rm temp;
