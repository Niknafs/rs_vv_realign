#!/usr/bin/perl
##*********************************************************************
##  rs_realignment_pipeline.v1.0.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-03-26
##*********************************************************************
## This script performs processing of bam alignments files including
#  .sorting 
#  .converting to fastq and formatting
#  .splitting for size
#  .re-alignment of fastq files with ELAND stand alone
#   
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_d $opt_o/;
getopts("d:o:");
my $usage =
".USAGE.   
rs_realignment_pipeline.v1.0.pl -d < dir of bam/bai files > -o < output directory >

.DESCRIPTION.
Realignment pipeline from bam to fastq back to bam

.OPTIONS.
  -d dir of bam/bai files
  -o output directory

.KEYWORDS.
eland, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_o;

my $BAM_DIR = $opt_d;
my $OUT_DIR = $opt_o;

# STAGE 0 -- setup
if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
  `mkdir $OUT_DIR/status`;
}

my %BAM_FILES = ();
my $BAM_LIST  = `ls sample_bams/*bam`;
chomp($BAM_LIST);
my @BAM_LIST  = split "\n", $BAM_LIST;
foreach my $bam (@BAM_LIST){
  my @bam = split /\//, $bam;
  $BAM_FILES{$bam[$#bam]} = $bam;
}

foreach my $bam (sort keys %BAM_FILES){
  print "$bam\n";
  # samtools sort sub.bam sub.sorted
  # bedtools bamtofastq -i sub.sorted.bam -fq sub.sorted.R1.fastq -fq2 sub.sorted.R2.fastq
  # split -l 20000 -d sub.sorted.R1.fastq sub.sorted.R1.fastq.
  # split -l 20000 -d sub.sorted.R2.fastq sub.sorted.R2.fastq.
  #
  # sub.sorted.R1.fastq.00
  # sub.sorted.R1.fastq.01
  # sub.sorted.R1.fastq.02....
  #
  # sed 's/\// /g;n;n;n;' sub.sorted.R1.fastq
  # sed 's/\// /g;n;n;n;' sub.sorted.R2.fastq

  
}

