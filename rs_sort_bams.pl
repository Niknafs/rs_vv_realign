#!/usr/bin/perl
##*********************************************************************
##  rs_sort_bams.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-02
##*********************************************************************
## This script performs processing of bam alignments files including
#  .sorting 
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_d $opt_o $opt_s/;
getopts("d:o:s:");
my $usage =
".USAGE.   
rs_sort_bams.pl -d < dir of bam files > -o < output directory >

.DESCRIPTION.

.OPTIONS.
  -d   dir of input bam files
  -o   output directory

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
  `chmod ug+rw $OUT_DIR`;
}

if (! -e "$OUT_DIR/status"){
  `mkdir $OUT_DIR/status`;
  `chmod ug+rw $OUT_DIR/status`;
}

my %BAM_FILES = ();
my $BAM_LIST  = `ls $BAM_DIR/*bam`;
chomp($BAM_LIST);
my @BAM_LIST  = split "\n", $BAM_LIST;
foreach my $bam (@BAM_LIST){
  my @bam = split /\//, $bam;
  $BAM_FILES{$bam[$#bam]} = $bam;
}

foreach my $bam (sort keys %BAM_FILES){
  print "$bam...\n";
  my $fp_bam = $BAM_FILES{$bam};
  my $prefix = $bam;
  $prefix =~ s/\.bam//g;
   
  my $ANALYSIS_DIR = "$OUT_DIR/analysis\_$prefix";

  if (-e "$OUT_DIR/status/$prefix.bam_to_fq.ck"){
    next;
  }else{
    `echo processing > $OUT_DIR/status/$prefix.bam_to_fq.ck`;
    `mkdir $ANALYSIS_DIR`; 
  }

  # SORT BAM FILE
  # begin by screening out poor quality sequence
  `samtools sort $fp_bam $ANALYSIS_DIR/$prefix.sorted`;
  `samtools index $ANALYSIS_DIR/$prefix.sorted.bam`; 
}
