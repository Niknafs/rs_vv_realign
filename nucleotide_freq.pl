#!/usr/bin/perl
#*********************************************************************
#  rs_coverage.pl*
#  Author:  James Robert White, PhD
#  Email:   james.dna.white@gmail.com
#  Created: 2015-10-22
#*********************************************************************

#*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use File::Spec;
use warnings;
use strict;
use Switch;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/;
#*********************************************************************

use vars qw/$opt_o $opt_b $opt_i/;
getopts("i:b:o:");
my $usage =
".USAGE.
rs_bam_to_nucleotidefreqs.pl -i < sorted indexed bam file > -b < input bed file of regions > -o < output directory >

TO RUN INTERACTIVELY YOU MUST REQUEST A LARGE MEM_FREE RESOURCE E.G.
qrsh -l mf=40G,h_vmem=40G -l gwas -q gwas.q\@compute-059 -pe local 1

NOTE THE INPUT BAM FILE MUST BE SORTED AND INDEXED

NOTE THE FOLLOWING MODULES ARE REQUIRED
module load perl
module load samtools/1.1
module load bedtools

.KEYWORDS.
coverage calculation
\n";

die $usage unless defined $opt_i
              and defined $opt_b
              and defined $opt_o;

my $bam     = File::Spec->rel2abs( $opt_i );
my $bed     = $opt_b;
my $outDir  = $opt_o;

print "Input bam: $bam\n";
print "Input bam: $bed\n";
print "Specified output directory: $outDir\n";

# PROCESSING ****************************************************************
# PARAMETERS
my $PICARD_DIR = "~jrwhite/jlib/picard-tools-1.84";

if (! -e "$outDir"){
  `mkdir $outDir`;
  `chmod ug+rw $outDir`;
}

# else{
#  die "Directory $outDir already exists...\n";
# }

# Begin by reheadering the bam artificially increasing the header lines
if (! -e "$outDir/header.sam"){
    `samtools view -H $bam > $outDir/header.sam`;
    open IN, "$outDir/header.sam" or die;
    open OUT, ">$outDir/reheader.sam" or die;
    while(<IN>){
        chomp($_);
        if ($_ =~ /^\@SQ/){
            my @A = split /\tLN\:/, $_;
            $A[1] = $A[1]+1000;
            $_ = join("\tLN:", @A);
        }
        print OUT "$_\n";
    }
    close IN;
    close OUT;

    # REHEADER
    `samtools reheader $outDir/reheader.sam $bam > $outDir/orig_reheader.bam`;
}

if (! -e "$outDir/cleaned.bam" ){
    # Clean up a messy bam file
    `java -Xms8g -Xmx32g -XX:ParallelGCThreads=1 -jar $PICARD_DIR/CleanSam.jar I=$outDir/orig_reheader.bam O=$outDir/cleaned.bam`;
}

# Mark Duplicates and remove
if (! -e "$outDir/collapsed.bam" ){
    `java -Xms8g -Xmx32g -jar $PICARD_DIR/MarkDuplicates.jar I=$outDir/cleaned.bam O=$outDir/collapsed.bam METRICS_FILE=$outDir/collapsed.bam.metrics ASSUME_SORTED=True REMOVE_DUPLICATES=True CREATE_INDEX=True`;
}

# To extracted mapped reads where both mates mapped:
`samtools view -b -F12 $outDir/collapsed.bam > $outDir/collapsed.map-map.bam`;

# To extract mapped reads whose mates are unmapped:
`samtools view -b -F4 -f8 $outDir/collapsed.bam > $outDir/collapsed.map-unmap.bam`;

# MPILEUP
open MM, "samtools mpileup -l $bed -d 1e6 -AB $outDir/collapsed.map-map.bam |";
open MMO, ">$outDir/map-map.coverage.txt";
print MMO "Chromosome\tCoord\tTotal\tA\tC\tG\tT\n";
while (<MM>){
  # chr10.fa        54860   N       2       gG      i~
  my @A = split "\t", $_;
  my $str = uc($A[4]);
  my @str = split "", $str;
  my %counts = (
    "A" => 0,
    "C" => 0,
    "G" => 0,
    "T" => 0
  );
  foreach my $s (@str){
    $counts{$s}++;
  }
  print MMO "$A[0]\t$A[1]\t$A[3]\t$counts{A}\t$counts{C}\t$counts{G}\t$counts{T}\n";  
}  
close MM;
close MMO;

open MU, "samtools mpileup -l $bed -d 1e6 -AB $outDir/collapsed.map-unmap.bam |";
open MUO, ">$outDir/map-unmap.coverage.txt";
print MUO "Chromosome\tCoord\tTotal\tA\tC\tG\tT\n";
while (<MU>){
  # chr10.fa        54860   N       2       gG      i~
  my @A = split "\t", $_;
  my $str = uc($A[4]);
  my @str = split "", $str;
  my %counts = (
    "A" => 0,
    "C" => 0,
    "G" => 0,
    "T" => 0
  );
  foreach my $s (@str){
    $counts{$s}++;
  }
  print MUO "$A[0]\t$A[1]\t$A[3]\t$counts{A}\t$counts{C}\t$counts{G}\t$counts{T}\n";
}
close MU;
close MUO;

# END
