#!/usr/bin/perl
##*********************************************************************
##  rs_bam_to_fq.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-02
##*********************************************************************
## This script performs processing of bam alignments files including
#  .sorting 
#  .converting to fastq and formatting
#  .splitting for size
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
use vars qw/$opt_d $opt_o $opt_s/;
getopts("d:o:s:");
my $usage =
".USAGE.   
rs_bam_to_fq.pl -d < dir of bam/bai files > -o < output directory > -s < split final fastq size >

.DESCRIPTION.

.OPTIONS.
  -d   dir of bam/bai files
  -o   output directory
  -s   size to split the final fqs into (default 1e6)

.KEYWORDS.
eland, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_o;

my $BAM_DIR = $opt_d;
my $OUT_DIR = $opt_o;
my $ENTRIES_PER_ALIGNMENT = 1000000;
if (defined($opt_s)){
  $ENTRIES_PER_ALIGNMENT = $opt_s;
}
my $FASTQ_SPLIT_DEX = $ENTRIES_PER_ALIGNMENT*4;

# STAGE 0 -- setup
if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
}

if (! -e "$OUT_DIR/status"){
  `mkdir $OUT_DIR/status`;
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
  `samtools sort -n $fp_bam $ANALYSIS_DIR/$prefix.sorted`;

  # CONVERT TO PAIRED END FASTQ
  `bedtools bamtofastq -i $ANALYSIS_DIR/$prefix.sorted.bam -fq $ANALYSIS_DIR/tmp.R1.fastq -fq2 $ANALYSIS_DIR/tmp.R2.fastq`;
  `sed 's/\\// /g;n;n;n;' $ANALYSIS_DIR/tmp.R1.fastq > $ANALYSIS_DIR/$prefix.R1.fastq`;
  `sed 's/\\// /g;n;n;n;' $ANALYSIS_DIR/tmp.R2.fastq > $ANALYSIS_DIR/$prefix.R2.fastq`;
  
  # clean up
  `rm $ANALYSIS_DIR/tmp.R1.fastq`;
  `rm $ANALYSIS_DIR/tmp.R2.fastq`;

  # clean up
  `rm $ANALYSIS_DIR/$prefix.sorted.bam`;
  
  `split -l $FASTQ_SPLIT_DEX -a 3 -d $ANALYSIS_DIR/$prefix.R1.fastq $ANALYSIS_DIR/$prefix.R1.fastq`;
  `split -l $FASTQ_SPLIT_DEX -a 3 -d $ANALYSIS_DIR/$prefix.R2.fastq $ANALYSIS_DIR/$prefix.R2.fastq`;
  
  # clean up 
  `rm $ANALYSIS_DIR/$prefix.R1.fastq`;
  `rm $ANALYSIS_DIR/$prefix.R2.fastq`;

  # how many indices were generated?
  my $r1_n = `ls $ANALYSIS_DIR/$prefix.R1.fastq* | wc | awk {'print \$1'}`;
  chomp($r1_n);
  $r1_n--;
  
  for my $i (0 .. ($r1_n-1)){
    my $dex = sprintf("%03d", $i);
    # $ANALYSIS_DIR/$prefix.R1.fastq$dex
    # $ANALYSIS_DIR/$prefix.R2.fastq$dex
 
    #print "ELAND FOR DEX $dex...\n"; 
    # `ELAND_standalone.pl -if $ANALYSIS_DIR/$prefix.R1.fastq$dex -if $ANALYSIS_DIR/$prefix.R2.fastq$dex --use-bases y100 --use-bases y100 -ref /tmp/jrwhite -it FASTQ -od $ANALYSIS_DIR/eland_out\_$dex --bam`; 
  } 
}