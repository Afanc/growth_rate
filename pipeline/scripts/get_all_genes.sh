#!/bin/bash                                                                                                            

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="get all genes"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

for file in `ls gff/*.gff`; do
    id=`echo $file | cut -d. -f1 | cut -d"/" -f2`;
    #echo "working on "$file", output is core/"$id".core";
    cat $file | awk 'BEGIN{OFS=" "} $3=="gene"{print $4,$5,$7,$9}' | grep gene= | awk 'BEGIN{FS=";"; OFS=" "} {print $1, $3}' | awk 'BEGIN{FS=" ";}{print $5,$1,$2,$3}' | sed 's/gene=//' > "full_genome/"$id".gen";
done
