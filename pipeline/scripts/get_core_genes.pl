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

my $empty = 1;

my $out_file = 'non_core_genes.txt';
open(my $out, '>', $out_file) or die "Could not open file '$out_file' $!";

my @gff = glob("gff/*.gff");
#(!@gff) ? die "no gff files present !" : print "all good, got gff files\n";
die "no gff files present !" if (!@gff);
open (genes_list, "joined_core") or die "can't read joined_core" ;

while (my $gene = <genes_list>) {
    chomp $gene;
    #print("working on ", $gene, "\n");

    foreach my $file (@gff) {
        open(my $f, "<", $file) or die "can't read file $!";

        my $exists_in_genome = 0;

        while (<$f>){
            if ($_ =~ m/$gene/ig){
                $exists_in_genome = 1;
                last;
            }
        }
        close($f);
        if (!$exists_in_genome){
            print $out "\t> gene ", $gene, " is not in all genomes ! \n";
            $empty=0;
            last;
        }
    }
}

if ($empty){
    print $out "all genes present in all strains !\n";
}
close(genes_list);
close($out);
