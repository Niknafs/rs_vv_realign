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
#
#*********************************************************************
#
# Notes regarding SGE settings:
#
#   - Both the merging and sorting operations can use multiple cores (e.g.,
#     specify 'pe local 16'). The overhead for RAM is fairly small for merging.
#     Thousands of ELAND bam files from an ELAND realignment can be merged using
#     less than 5G.
#
#   - The sorting can require a very large amount of RAM. By default, the
#     bamfiles are sorted with a maximum allowed memory of 10G per thread. The
#     following SGE settings would require a total of approx. 150G RAM:
#     -pe local 16
#     -l mem_free=15G
#     -l h_vmem=16G
#
#   - If the sort fails, delete the relevant status/ files and submit the SGE
#     script with the above settings. This script will recognize that the merged
#     file has already been created and proceed directly to the sorting.
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
use vars qw/$opt_d $opt_i $opt_o $opt_s/;
getopts("d:i:o:");
my $usage =
".USAGE.   
rs_merge_sample_bams.pl -d < dir w/ reanalysis.bam files > -i < file index number to process > -o < final output directory > 

.OPTIONS.
  -d  post alignment ouput directory of files
      Here we assume there exist subdirectories with sample_ids 
      as prefixes (e.g. <prefix>.XXX_<eland or other lowercase suffix>)
  -i  file index number to process
  -o  final output directory (e.g. merged_bams)
\n";

die $usage unless defined $opt_d
              and defined $opt_i
              and defined $opt_o;

my $F_DIR   = $opt_d;
my $INDEX   = $opt_i;
my $OUT_DIR = $opt_o;

# STAGE 0 -- setup
if (! -e $F_DIR){
  die "Error $F_DIR not found!\n";
}

if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
  `mkdir $OUT_DIR/status`;
  `chmod ug+rw $OUT_DIR`;
  `chmod ug+rw $OUT_DIR/status`;  
}

if (! -e "$OUT_DIR/status"){
  `mkdir $OUT_DIR/status`;  
  `chmod ug+rw $OUT_DIR/status`;
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
  $s[$#s] =~ s/\.[0-9][0-9][0-9][0-9]\_[a-z]+$//g;
  $SAMPLE_LIST{$s[$#s]} = 1;
}

my @id_list = sort keys %SAMPLE_LIST;
if (($INDEX-1) > $#id_list){
  die "No corresponding file index for $INDEX\n";
}

my $id      = $id_list[($INDEX-1)];
if (-e "$OUT_DIR/status/$id.merge.ck"){
  print "$id has already been processed...\n";
  exit;
}else{
  print "Processing $id...\n";
  `echo processing > $OUT_DIR/status/$id.merge.ck`;
}

# merge those bams
print "Begin merge for $id...\n";
`samtools merge -@ 15 $OUT_DIR/$id.merged.bam \`find $F_DIR/$id\.* -name "reanalysis.bam" -print\``;
print "Begin sort for $id...\n";
`samtools sort -m 10G -@ 15 $OUT_DIR/$id.merged.bam $OUT_DIR/$id.sorted`;
print "Begin index for $id...\n";
`samtools index $OUT_DIR/$id.sorted.bam`;
# `mv $OUT_DIR/$id.sorted.bam.bai $OUT_DIR/$id.sorted.bai`; 

# print "Removing original merged bam for $id...\n";
# `rm $OUT_DIR/$id.merged.bam`;
