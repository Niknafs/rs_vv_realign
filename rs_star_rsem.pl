#!/usr/bin/perl
##*********************************************************************
##  rs_transcriptomealignment.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2016-04-28
##*********************************************************************
## This script performs processing of fastq.gz files
## followed by alignment using HiSat2
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_d $opt_o $opt_i $opt_r $opt_c/;
getopts("d:r:i:o:c:");
my $usage =
".USAGE.   
rs_star_rsem.pl -d < paired end fastq.gz files > -i < file index number to process > -r < path to reference directory > -o < output dir >

.DESCRIPTION.

.OPTIONS.
  -d  paired end fastq.gz files
  -i  index number of the file to process
  -r  path to reference transcriptome index (prefix of index files)
  -c  number of cores to use (optional - default 16)
  -o  output directory for alignment files
  
.KEYWORDS.
star, rsem, transcriptome, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_i
              and defined $opt_r
              and defined $opt_o;

my $FQ_DIR  = $opt_d;
my $INDEX   = $opt_i;
my $REF_DIR = $opt_r;
my $OUT_DIR = $opt_o;
my $max_cores = 16;

if (defined $opt_c){
  $max_cores = $opt_c;
}

# STAGE 0 -- setup
if (! -e $FQ_DIR){
  die "Error $FQ_DIR not found!\n";
}

if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
  `chmod ug+rw $OUT_DIR`;
}
if (! -e "$OUT_DIR/status"){
  `mkdir $OUT_DIR/status`;
  `chmod ug+rw $OUT_DIR/status`;
}
if (! -e "$OUT_DIR/rsem"){
  `mkdir $OUT_DIR/rsem`;
  `chmod ug+rw $OUT_DIR/rsem`;
}

my %FQ_FILES = ();
my $FQ_LIST  = `find $FQ_DIR -name "*.fastq.gz" -print`;
chomp($FQ_LIST);
my @FQ_LIST  = split "\n", $FQ_LIST;
my %FQ_PAIRS = ();

# create the paired end list on the fly
foreach my $fq (@FQ_LIST){
  my @fq     = split /\//, $fq;
  my $id = $fq[$#fq];
  $id    =~ s/\_R[12]\.fastq\.gz//g;
  my $dir = join("/", @fq[0 .. ($#fq-1)]);

  $FQ_PAIRS{$id}{"DIR"} = $dir;
  if ($fq[$#fq] =~ /\_R2.fastq.gz/){
    $FQ_PAIRS{$id}{"R2"} = $fq;  
  }else{ 
    $FQ_PAIRS{$id}{"R1"} = $fq; 
  } 
}

my @id_list = sort keys %FQ_PAIRS;

if (($INDEX-1) > $#id_list){
  die "No corresponding file index for $INDEX\n";
}

my $id      = $id_list[($INDEX-1)];
if (-e "$OUT_DIR/status/$id.rsem.ck"){
  die "$id has already been processed according to its status file.\n";
}else{
  `echo processing > $OUT_DIR/status/$id.rsem.ck`;
}

if (!defined($FQ_PAIRS{$id}{"R1"}) or !defined($FQ_PAIRS{$id}{"R2"})){
  die "Missing R1 or R2 pair for $id\n";
}
 
print "STAR-RSEM RUNS FOR $id ...\n"; 
`mkdir $OUT_DIR/rsem/$id`;

print "rsem-calculate-expression --paired-end --star --star-path /users/jrwhite/jbin --star-gzipped-read-file --sort-bam-by-coordinate --output-genome-bam -p $max_cores $FQ_PAIRS{$id}{R1} $FQ_PAIRS{$id}{R2} $REF_DIR $OUT_DIR/rsem/$id/$id";
 `rsem-calculate-expression --paired-end --star --star-path /users/jrwhite/jbin --star-gzipped-read-file --sort-bam-by-coordinate --output-genome-bam -p $max_cores $FQ_PAIRS{$id}{R1} $FQ_PAIRS{$id}{R2} $REF_DIR $OUT_DIR/rsem/$id/$id`;
