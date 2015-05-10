#!/usr/bin/perl
##*********************************************************************
##  rs_eland.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-02
##*********************************************************************
## This script performs processing of split fastq files including
#  with ELAND to generate aligned bam files
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
getopts("d:r:");
my $usage =
".USAGE.   
rs_eland.pl -d < dir of bam_to_fq output files > -r < path to reference genome fasta >

.DESCRIPTION.

.OPTIONS.
  -d  post bam_to_fq ouput directory of files (subdirectories all start with analysis)
  -r  path to reference genome fasta file 
  

.KEYWORDS.
eland, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_r;

my $FQ_DIR  = $opt_d;
my $REF_DIR = $opt_r;

# STAGE 0 -- setup
if (! -e $FQ_DIR){
  die "Error $FQ_DIR not found!\n";
}
if (! -e $REF_DIR){
  die "Error $REF_DIR not found!\n";
}

my %FQ_FILES = ();
my $FQ_LIST  = `find $FQ_DIR -name "*.fastq*" -print`;
chomp($FQ_LIST);
my @FQ_LIST  = split "\n", $FQ_LIST;
my %FQ_PAIRS = ();

foreach my $fq (@FQ_LIST){
  next if ($fq =~ /fastq$/ or $fq =~ /fastq.gz$/);
  my @fq     = split /\//, $fq;
  my $id = $fq[$#fq];
  $id    =~ s/\.R[12]\.fastq/\./g;
  my $dir = join("/", @fq[0 .. ($#fq-1)]);

  $FQ_PAIRS{$id}{"DIR"} = $dir;
  if ($fq[$#fq] =~ /\.R2.fastq/){
    $FQ_PAIRS{$id}{"R2"} = $fq;  
  }else{ 
    $FQ_PAIRS{$id}{"R1"} = $fq; 
  } 
}

foreach my $id (sort keys %FQ_PAIRS){
   
  if (-e "$FQ_DIR/status/$id.ck"){
    next;
  }else{
    `echo processing > $FQ_DIR/status/$id.ck`;
  }

  if (!defined($FQ_PAIRS{$id}{"R1"}) or !defined($FQ_PAIRS{$id}{"R2"})){
   next;
  }
 
  print "ELAND FOR $id ...\n"; 
  `ELAND_standalone.pl -if $FQ_PAIRS{$id}{R1} -if $FQ_PAIRS{$id}{R2} --use-bases y100 --use-bases y100 -ref $REF_DIR -it FASTQ -od $FQ_PAIRS{$id}{DIR}/$id\_eland --bam`; 

  # clean up
  `rm $FQ_PAIRS{$id}{DIR}/$id\_eland/*.gz`;
  `rm $FQ_PAIRS{$id}{DIR}/$id\_eland/*.xml`;
  `rm $FQ_PAIRS{$id}{DIR}/$id\_eland/*.oa`;
  `rm $FQ_PAIRS{$id}{DIR}/$id\_eland/*.txt`;  
}
