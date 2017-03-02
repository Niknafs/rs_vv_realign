#!/usr/bin/perl
##*********************************************************************
##  rs_merge_sample_exports.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-10-15
##*********************************************************************
## This script performs 
# - identification and merging of partitioned bam file
# that are realigned
# - merging export files from ELAND runs
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
rs_merge_sample_exports.pl -d < dir w/ reanalysis_NoIndex_L001_R{1/2}_001_export.txt.gz files > -i < file index number to process > -o < final output directory > 

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

# merge those export.txt.gz R1/R2 files
print "Begin export merge for $id R1...\n";
`cat \`find $F_DIR/$id\.* -name "reanalysis_NoIndex_L001_R1_001_export.txt.gz" -print\` > $OUT_DIR/$id\_R1_001_export.txt.gz`;
print "Begin export merge for $id R2...\n";
`cat \`find $F_DIR/$id\.* -name "reanalysis_NoIndex_L001_R2_001_export.txt.gz" -print\` > $OUT_DIR/$id\_R2_001_export.txt.gz`;

