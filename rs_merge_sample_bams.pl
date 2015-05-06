#!/usr/bin/perl
##*********************************************************************
##  rs_merge_sample_bams.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-16
##*********************************************************************
## This script performs 
# - identification and merging of partitioned bam file
# that are realigned using ELAND
# - these final bams are deposited into a different (more permanent)
# directory
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
getopts("d:o:");
my $usage =
".USAGE.   
rs_merge_sample_bams.pl -d < dir of bam_to_fq output files > -o < final output directory > 

.DESCRIPTION.

.OPTIONS.
  -d  post bam_to_fq ouput directory of files (subdirectories all start with analysis)
  -o final output directory 

.KEYWORDS.
eland, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_o;

my $FQ_DIR  = $opt_d;
my $OUT_DIR = $opt_o;

# STAGE 0 -- setup
if (! -e $FQ_DIR){
  die "Error $FQ_DIR not found!\n";
}

if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
}

# *****************************************************
# underneath the FQ directory are analysis_*
# subdirectories
# each sample name needs to be extracted here
# *****************************************************
my $SAMPLE_LIST = `ls -d $FQ_DIR/analysis_*`;
chomp($SAMPLE_LIST);
my @SAMPLE_LIST = split "\n", $SAMPLE_LIST; 

my %SAMPLE_LIST = ();
foreach my $s (@SAMPLE_LIST){
  my @s = split /\//, $s;
  $s[$#s] =~ s/analysis_//g;
  $SAMPLE_LIST{$s[$#s]} = 1;
}

foreach my $id (sort keys %SAMPLE_LIST){
   
  if (-e "$FQ_DIR/status/$id.merge.ck"){
    next;
  }else{
    print "Processing $id...\n";
    `echo processing > $FQ_DIR/status/$id.merge.ck`;
  }

  # merge those bams
  # my $list = `find $FQ_DIR/analysis_$id -name "reanalysis.bam" -print`;
  `samtools merge $OUT_DIR/$id.merged.bam \`find $FQ_DIR/analysis_$id -name "reanalysis.bam" -print\``; 
}
