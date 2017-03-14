#!/usr/bin/perl
##*********************************************************************
##  rs_process_bams.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
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
use vars qw/$opt_d $opt_o $opt_i/;
getopts("d:i:o:");
my $usage =
".USAGE.   
rs_process_bams.pl -d < dir of bam files > -i < file index number to process > -o < output directory >

.DESCRIPTION.

.OPTIONS.
  -d   dir of input bam files
       -- currently all files in this directory ending in \".bam\" will be processed
       -- these bam files may contain \".fa\" naming scheme which 
       -- this script will correct those .fas and also sort and index
          the downstream bam files, producing <prefix>.sorted.bam 
  -i   index number of the file to process
       -- in a qsub script, the -t parameter will be used to supply which 
          index to use. e.g. -t 1-10 will be an array job and process bam
          files 1 through 10 (alphabetically)
  -o   output directory

.EXAMPLE QUB SCRIPT.
#!/bin/bash
#\$ -cwd
#\$ -j y
#\$ -l mem_free=2.5G
#\$ -l h_vmem=20G
#\$ -l h_rt=60:00:00
#\$ -pe local 16
#\$ -t 1-6
module load perl
rs_process_bam.pl -d \$PWD/data -i \${SGE_TASK_ID} -o \$PWD/processed

\n";

die $usage unless defined $opt_d
              and defined $opt_i
              and defined $opt_o;

my $BAM_DIR     = $opt_d;
my $INDEX       = $opt_i;
my $OUT_DIR     = $opt_o;
my $SAMTOOLS    = "/users/jrwhite/jlib/samtools-1.3/samtools";
my $max_threads = 16;

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
my $BAM_LIST  = `ls $BAM_DIR/*.bam`;
chomp($BAM_LIST);
my @BAM_LIST  = split "\n", $BAM_LIST;
foreach my $bam (@BAM_LIST){
  my @bam = split /\//, $bam;
  $BAM_FILES{$bam[$#bam]} = $bam;
}

my @id_list = sort keys %BAM_FILES;

if (($INDEX-1) > $#id_list){
  die "No corresponding file index for $INDEX\n";
}

my $id      = $id_list[($INDEX-1)];
if (-e "$OUT_DIR/status/$id.process.ck"){
  die "$id has already been processed according to its status file.\n";
}else{
  `echo processing > $OUT_DIR/status/$id.process.ck`;
}

print "$id...\n";
`echo removing fas >> $OUT_DIR/status/$id.process.ck`;
my $fp_bam = $BAM_FILES{$id};
my $prefix = $id;
$prefix =~ s/\.bam//g;
 
# generate nofas version
my $nofabam = "$OUT_DIR/$prefix.nofa.bam"; 
`$SAMTOOLS view -H $fp_bam | sed -e 's/.fa//' | $SAMTOOLS reheader - $fp_bam > $nofabam`; 

# SORT BAM FILE
# begin by screening out poor quality sequence
`echo sorting >> $OUT_DIR/status/$id.process.ck`;
`$SAMTOOLS sort -\@ $max_threads -T $OUT_DIR/$prefix.tmp $nofabam -o $OUT_DIR/$prefix.sorted.bam`;

`echo indexing >> $OUT_DIR/status/$id.process.ck`;
`$SAMTOOLS index $OUT_DIR/$prefix.sorted.bam`; 

# clean up
`rm $nofabam`; 

`echo complete. >> $OUT_DIR/status/$id.process.ck`;
