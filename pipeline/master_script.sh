#!/bin/bash                                                                                                            

#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="master"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

#rm bad.error good.out;
#touch bad.error;

echo -e "let's go !\n">good.out;
# foolproofing
if [[ -f repl_rates.csv && -f joined_core ]];
    then echo -e "all files are here\n">>good.out;
    else
        echo -e "either repl_rates.csv or joined_core is missing\n">>bad.error;
        exit 126;
fi

if [[ `head joined_core` == gene=* ]];
    then echo -e "joined_core seems ok\n">>good.out; 
    else cp .joined_core.bak joined_core;
fi

# arguments and stuff
PROGNAME=$0

usage() {
    cat << EOF >&2
Usage: $PROGNAME [-f <features>] [-o <output>] 

-p <params>: comma separated, no space ! [rank,beg,end,dir,occ]
-f <file>: name of .zip file [default is results.zip] 
-- anything else : will show this or will be ignored !
-- also, if you get the features wrong, default is rank,dir,occ
EOF
echo -e "\n>>>> SORRY DUDE, THERE WERE ERRORS\nCheck bad.error\n">>good.out;
exit 1
}

output_zip=results.zip features=rank,dir,occ
#while getopts f:o:v o; do
while getopts 'f:p:' o; do
    case $o in
        (p) features=$OPTARG;;
        (f) output_zip=$OPTARG;;
        (*) usage
    esac
done
shift "$((OPTIND - 1))"

features=`echo $features | tr ',' ' '`

echo -e "output file is "$output_zip" and paramaters are "$features"\n">>good.out;

# scripts
tot=7

echo -e "step 1/"$tot", get non_core genes\n">>good.out;
./scripts/get_core_genes.pl >> bad.error;

echo -e "step 2/"$tot", filter out non_core genes\n">>good.out;
./scripts/clean_core_genes.sh >> bad.error;

echo -e "step 3/"$tot", create core files\n">>good.out;
mkdir -p core; 
./scripts/run_filter_core_genes.sh >> bad.error;

echo -e "step 4/"$tot", create stats files\n">>good.out;
mkdir -p stats;
./scripts/extract_infos.pl >> bad.error;

echo -e "step 5/"$tot", pca and kmeans\n">>good.out;
mkdir -p results;
/home/dmollet/.conda/envs/python_env/bin/python3 ./scripts/pca_kmeans.py $features >> bad.error;

echo -e "step 6/"$tot", correlations\n">>good.out;
/home/dmollet/.conda/envs/python_env/bin/python3 ./scripts/correlations.py $features >> bad.error;

echo -e "step 7/"$tot", zipping everything and prepare for upload\n">>good.out;
rm -f "results/"$output_zip;
zip -r "results/"$output_zip results/*;
mv "results/"$output_zip ~/;

echo -e "done ! \n"$output_zip" is ready in your home\n">>good.out;

if [ -s bad.error ]
    then
        echo -e "\n>>>> SORRY DUDE, THERE WERE ERRORS\nCheck bad.error\n">>good.out;
fi
