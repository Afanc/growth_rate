#!/usr/bin/perl
#SBATCH --mail-user=dariush.mollet@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="filter"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=48G
#SBATCH --output=good.out
#SBATCH --error=bad.error

use warnings;
use strict; 

my @core = glob("full_genome/*.gen");

foreach my $file_in (@core) { #each core file at a time

    my $file_out = "stats_all/".substr($file_in, 12 , -4).".stats_all";
    #print($file_out,"\n");

    open(my $out, ">", $file_out) or die "can't write to file $!\n";
    open (genes_list, "all_genes.txt") or die "can't read all_genes.txt\n" ;
    print $out "gene rank beg end dir occ\n";

    while (my $gene = <genes_list>) { #open genes' list

        open(my $in, "<", $file_in) or die "can't read file $!\n";
        chomp $gene;
        my $occurances = 0;

        print $out $gene." ";
        my $first = 1;
        while (<$in>){  #go through core file for each gene and get occurences
            if ($_ =~ m/$gene/g){
                my @line = split(" ", $_);
                my $direction = ($line[3] eq '+') ? 1 : 0;      #1 if +, 0 else
                if ($first){
                    print $out $.,  " ", $line[1], " ", $line[2], " ", $direction, " ";
                    $first = 0;
                }
                $occurances += 1;

                #test for all strains being assembled in the right direction
                if ($gene eq "dnaA" and $direction == 0){
                    die ">>>>>> OH SHIT, dude check ".$file_in." - it's heading in the wrong direction !\n";
                }
            }
        }

        if ($occurances != 0){
            print $out $occurances,"\n";
        }
        else {
            print $out 0, " ", 0, " ", 0, " ", 0, " ", 0, "\n";
        }
        close($in);
    }
    close(genes_list);
    close($file_out)
}
