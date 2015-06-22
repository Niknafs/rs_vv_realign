#!/usr/bin/perl
##*********************************************************************
##  rs_merge_sample_bams.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-16
##*********************************************************************
## This script performs 
# - identification and merging of partitioned bam file
# that are realigned
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
rs_merge_sample_bams.pl -d < dir w/ reanalysis.bam files > -o < final output directory > 

.OPTIONS.
  -d  post alignment ouput directory of files
      Here we assume there exist subdirectories with sample_ids 
      as prefixes (e.g. <prefix>.XXX_<eland or other lowercase suffix>)
  -o  final output directory (e.g. merged_bams)
\n";

die $usage unless defined $opt_d
              and defined $opt_o;

my $F_DIR  = $opt_d;
my $OUT_DIR = $opt_o;

# STAGE 0 -- setup
if (! -e $F_DIR){
  die "Error $F_DIR not found!\n";
}

if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
  `mkdir $OUT_DIR/status`;
}

if (! -e "$OUT_DIR/status"){
  `mkdir $OUT_DIR/status`;  
}

# *****************************************************
# underneath the F directory are <prefix>.XXX_<eland or other suffix>
# subdirectories
# each sample name needs to be extracted here
# *****************************************************
my $SAMPLE_LIST = `ls -d $F_DIR/*`;
chomp($SAMPLE_LIST);
my @SAMPLE_LIST = split "\n", $SAMPLE_LIST; 

my %SAMPLE_LIST = ();
foreach my $s (@SAMPLE_LIST){
  my @s = split /\//, $s;
  next if ($s[$#s] =~ /^status$/);
  $s[$#s] =~ s/\.[0-9][0-9][0-9]\_[a-z]+$//g;
  $SAMPLE_LIST{$s[$#s]} = 1;
}

foreach my $id (sort keys %SAMPLE_LIST){
   
  if (-e "$OUT_DIR/status/$id.merge.ck"){
    next;
  }else{
    print "Processing $id...\n";
    `echo processing > $OUT_DIR/status/$id.merge.ck`;
  }

  # merge those bams
  # my $list = `find $F_DIR/analysis_$id -name "reanalysis.bam" -print`;
  `samtools merge $OUT_DIR/$id.merged.bam \`find $F_DIR/$id\.* -name "reanalysis.bam" -print\``; 
}
