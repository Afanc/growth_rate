#!/bin/bash                                                                                                            

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="filter core genes"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

for file in `ls gff/*.gff`; do
    id=`echo $file | cut -d. -f1 | cut -d"/" -f2`;
    #echo "working on "$file", output is core/"$id".core";
    grep -iFf joined_core $file | awk 'BEGIN{FS=";"; OFS="\t"} FNR%2 {if ($4 ~/gene=|Name=|product=/) print $1, $4; else print $1, $5}' | awk 'BEGIN{FS="\t"}{print $10, $4, $5, $7}' | sed 's/^.*=//' > "core/"$id".core";
done
sed -i -e 's/gene=//' joined_core 
