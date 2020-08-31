#!/bin/bash                                                                                                            

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="clustal-om"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

export PATH=/software/bin:/home/dmollet/mcl/src/alien/oxygen/src/:$PATH;
module use /software/module/
module add vital-it
module add SequenceAnalysis/MultipleSequenceAlignment/clustal-omega/1.2.4

clustalo -i all_in_one.fasta -t DNA -v --outfmt=phy -o clustal_out

echo -e "--------\nyeyeye">>good.out
