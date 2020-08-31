#!/usr/bin/env perl

my $fastaFile = $ARGV[0];   # first param (well oriented fasta file)
my $origin = $ARGV[1]; #second param numeric value giving the position of oriC (ATG of dnaA)

use lib "$ENV{HOME}/lib";
use Bio::SeqIO;
my $seqin = Bio::SeqIO->new(-file => "$fastaFile", '-format' => 'Fasta');
my $seqout = Bio::SeqIO->new(-fh => \*STDOUT, '-format' => 'Fasta');

while(my $seq = $seqin->next_seq){
   $seq2=$seq->subseq($origin,$seq->length());
   $seq1=$seq->subseq(1,$origin-1);
   $seq2=Bio::Seq->new(-seq => "$seq2"."$seq1", -display_id => $seq->display_id());
          $seqout->write_seq($seq2);
}

exit 0; 
