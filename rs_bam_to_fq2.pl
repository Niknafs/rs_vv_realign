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
use vars qw/$opt_d $opt_i $opt_o $opt_s/;
getopts("d:i:o:s:");
my $usage =
".USAGE.   
rs_bam_to_fq.pl -d < dir of bam/bai files > -i < file index number to process > -o < output directory > -s < split final fastq size >

.DESCRIPTION.

.OPTIONS.
  -d   dir of bam/bai files
  -i  index number of the file to process
  -o   output directory
  -s   size to split the final fqs into (default 1e6)

.KEYWORDS.
eland, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_o;

my $BAM_DIR = $opt_d;
my $OUT_DIR = $opt_o;
my $INDEX   = $opt_i;
my $ENTRIES_PER_ALIGNMENT = 1000000;
if (defined($opt_s)){
  $ENTRIES_PER_ALIGNMENT = $opt_s;
}
my $FASTQ_SPLIT_DEX = $ENTRIES_PER_ALIGNMENT*4;

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

my @id_list = sort keys %BAM_FILES;
if (($INDEX-1) > $#id_list){
  die "No corresponding file index for $INDEX\n";
}
my $bam     = $id_list[($INDEX-1)];

print "$bam...\n";
my $fp_bam = $BAM_FILES{$bam};
my $prefix = $bam;
$prefix =~ s/\.bam//g;
   
my $ANALYSIS_DIR = "$OUT_DIR/analysis\_$prefix";

if (-e "$OUT_DIR/status/$prefix.bam_to_fq.ck"){
  die "$bam has already been processed according to its status file.\n";
}elsif(-e $ANALYSIS_DIR){
  die "WARNING/ERROR: $prefix has no status file ($prefix.bam_to_fq.ck), but it's analysis directory $ANALYSIS_DIR exists! $prefix will not be processed as a result\n";  
}else{
  `echo processing > $OUT_DIR/status/$prefix.bam_to_fq.ck`;
  `mkdir $ANALYSIS_DIR`; 
  `chmod ug+rw $ANALYSIS_DIR`;
}

# SORT BAM FILE
# begin by screening out poor quality sequence
# `samtools sort -m 32G -n $fp_bam $ANALYSIS_DIR/$prefix.sorted`;

# CONVERT TO PAIRED END FASTQ
# `bedtools bamtofastq -i $ANALYSIS_DIR/$prefix.sorted.bam -fq $ANALYSIS_DIR/tmp.R1.fastq -fq2 $ANALYSIS_DIR/tmp.R2.fastq`;
`/users/jrwhite/jlib/samtools-1.3/samtools fastq -1 $ANALYSIS_DIR/tmp.R1.fastq -2 $ANALYSIS_DIR/tmp.R2.fastq $fp_bam`;
`cat $ANALYSIS_DIR/tmp.R1.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $ANALYSIS_DIR/$prefix.R1.fastq`;
`cat $ANALYSIS_DIR/tmp.R2.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > $ANALYSIS_DIR/$prefix.R2.fastq`;
 
# clean up tmp fastqs aftersplitting
#`rm $ANALYSIS_DIR/tmp.R1.fastq`;
#`rm $ANALYSIS_DIR/tmp.R2.fastq`;

# clean up sorted bam file
#`rm $ANALYSIS_DIR/$prefix.sorted.bam`;
  
`split -l $FASTQ_SPLIT_DEX -a 4 -d $ANALYSIS_DIR/$prefix.R1.fastq $ANALYSIS_DIR/$prefix.R1.fastq`;
`split -l $FASTQ_SPLIT_DEX -a 4 -d $ANALYSIS_DIR/$prefix.R2.fastq $ANALYSIS_DIR/$prefix.R2.fastq`;
  
# clean up full files fastqs 
#`rm $ANALYSIS_DIR/$prefix.R1.fastq`;
#`rm $ANALYSIS_DIR/$prefix.R2.fastq`;
