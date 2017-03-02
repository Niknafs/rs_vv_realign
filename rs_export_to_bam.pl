#!/usr/bin/perl
##*********************************************************************
##  rs_export_to_bam.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-05-08
##*********************************************************************
## This script performs processing of CASAVA export files
#  to a sorted bam file
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_d $opt_o $opt_r/;
getopts("d:o:r:");
my $usage =
".USAGE.   
rs_export_to_bam.pl -d < dir containing export files > -r < reference fasta > -o < output directory >

.DESCRIPTION.

.OPTIONS.
  -d  directory containing subdirectories with one or more export files
      
      Paired Export files are expected to end in *_R1_export.txt or *_R2_export.txt 
      Unpaired Export files are expected to end in *export.txt

  -r  reference fasta file e.g. /dcs01/oncbio/rscharpf/White/jdb/hg18/human_hg18.fa
  -o  output directory with intermediate and final bam files


.KEYWORDS.
\n";

die $usage unless defined $opt_d 
              and defined $opt_r
              and defined $opt_o;

my $EX_DIR  = $opt_d;
my $OUT_DIR = $opt_o;
my $REF_FNA = $opt_r;

# STAGE 0 -- setup
if (! -e $EX_DIR){
  die "Error $EX_DIR not found!\n";
}

if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
  `mkdir $OUT_DIR/status`;
}

my %EX_FILES = ();
my $EX_LIST  = `find $EX_DIR -name "*export.txt" -print`;
chomp($EX_LIST);
my @EX_LIST  = split "\n", $EX_LIST;
my %EX_PAIRS = ();

foreach my $ex (@EX_LIST){
  my @ex = split /\//, $ex;
  my $id = $ex[$#ex];
  if ($id =~ /\_R[12]\_export.txt/){
    $id    =~ s/\_R[12]\_export.txt//g;
  }else{ # it's not a paired set
    $id =~ s/[\.\-\_\,\;\:]export.txt//g;
  }

  my $dir = join("/", @ex[0 .. ($#ex-1)]);

  $EX_PAIRS{$id}{"DIR"} = $dir;
  if ($ex[$#ex] =~ /R2_export.txt/){
    $EX_PAIRS{$id}{"R2"} = $ex;  
  }else{ 
    $EX_PAIRS{$id}{"R1"} = $ex; 
  } 
}

foreach my $id (sort keys %EX_PAIRS){
   
  if (-e "$OUT_DIR/status/$id.ck"){
    next;
  }else{
    `echo processing > $OUT_DIR/status/$id.ck`;
  }
   
  `mkdir $OUT_DIR/analysis_$id`;

  my $numFiles = 2;
  if (!defined($EX_PAIRS{$id}{"R1"}) or !defined($EX_PAIRS{$id}{"R2"})){
    $numFiles = 1;  
  }
 
  print "EXPORT_TO_SAM_TO_BAM FOR $id ...\n";
  if ($numFiles == 2){
    `/dcs01/oncbio/rscharpf/White/CASAVA-1.8.2/bin/illumina_export2sam.pl --read1=$EX_PAIRS{$id}{R1} --read2=$EX_PAIRS{$id}{R2} | samtools view -S -h -T $REF_FNA -b - | samtools sort - $OUT_DIR/analysis_$id/$id.sorted`; 
  }else{
    `/dcs01/oncbio/rscharpf/White/CASAVA-1.8.2/bin/illumina_export2sam.pl --read1=$EX_PAIRS{$id}{R1} | samtools view -S -h -T $REF_FNA -b - | samtools sort - $OUT_DIR/analysis_$id/$id.sorted`;
  }
  `samtools index $OUT_DIR/analysis_$id/$id.sorted.bam`;

  # clean up
}
