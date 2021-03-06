#!/usr/bin/perl
##*********************************************************************
##  rs_novoalign.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-02
##*********************************************************************
## This script performs processing of split fastq files including
#  with NOVOALIGN to generate aligned bam files
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_d $opt_o $opt_i $opt_r/;
getopts("d:r:i:o:");
my $usage =
".USAGE.   
rs_novoalign.pl -d < dir of bam_to_fq output files > -i < file index number to process > -r < path to reference genome novoalign index file (.nix or .ndx for example) > -o < output dir >

.DESCRIPTION.

.OPTIONS.
  -d  post bam_to_fq ouput directory of files (subdirectories all start with analysis)
  -i  index number of the file to process
  -r  path to reference genome novoalign index file (.nix or .ndx for example)
  -o  output directory for alignment files
  

.KEYWORDS.
alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_i
              and defined $opt_r 
              and defined $opt_o;

my $FQ_DIR  = $opt_d;
my $INDEX   = $opt_i;
my $REF_DIR = $opt_r;
my $OUT_DIR = $opt_o;

# STAGE 0 -- setup
if (! -e $FQ_DIR){
  die "Error $FQ_DIR not found!\n";
}
if (! -e $REF_DIR){
  die "Error $REF_DIR not found!\n";
}

if (! -e $OUT_DIR){
  `mkdir $OUT_DIR`;
  `chmod ug+rw $OUT_DIR`;
}
if (! -e "$OUT_DIR/status"){
  `mkdir $OUT_DIR/status`;
  `chmod ug+rw $OUT_DIR/status`;
}

my %FQ_FILES = ();
my $FQ_LIST  = `find $FQ_DIR -name "*.fastq*" -print`;
chomp($FQ_LIST);
my @FQ_LIST  = split "\n", $FQ_LIST;
my %FQ_PAIRS = ();

foreach my $fq (@FQ_LIST){
  next if ($fq =~ /fastq$/ or $fq =~ /fastq.gz$/);
  my @fq  = split /\//, $fq;
  my $id  = $fq[$#fq];
  $id     =~ s/\.R[12]\.fastq/\./g;
  my $dir = join("/", @fq[0 .. ($#fq-1)]);

  $FQ_PAIRS{$id}{"DIR"} = $dir;
  if ($fq[$#fq] =~ /\.R2.fastq/){
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

if (-e "$OUT_DIR/status/$id.novoalign.ck"){
  die "$id has already been processed according to its status file.\n";
}else{
  `echo processing > $OUT_DIR/status/$id.novoalign.ck`;
}

if (!defined($FQ_PAIRS{$id}{"R1"}) or !defined($FQ_PAIRS{$id}{"R2"})){
   die "Missing R1 or R2 pair for $id\n";
}
 
print "NOVOALIGN FOR $id...\n";
if (-e "$OUT_DIR/$id\_novo"){
  die "$id has already been processed according to it's present directory with outputs\n";
}

`mkdir $OUT_DIR/$id\_novo`;
`(time novoalign -d $REF_DIR -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -f $FQ_PAIRS{$id}{R1} $FQ_PAIRS{$id}{R2} -oSAM -oFullNW | samtools view -bS - > $OUT_DIR/$id\_novo/reanalysis.bam) 2>$OUT_DIR/$id\_novo/time.txt`;

