#!/bin/bash                                                                                                            

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="trnas"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

module add UHTS/Analysis/aragorn/1.2.38

for file in `ls fna/*.fna`; do
    id=`echo $file | cut -d. -f1 | cut -d"/" -f2`;
    #echo "working on "$file", output is core/"$id".core";
    aragorn -l -gc4  -w $file | grep -v 'tmRNA' | awk '{print $5}' | tr -d '()' | tail -n+3 | sort | uniq -c | awk 'BEGIN{OFS=",";} {print $2, $1}' > "trna/"$id".trna"
done


