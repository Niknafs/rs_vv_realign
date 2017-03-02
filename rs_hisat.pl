#!/usr/bin/perl
##*********************************************************************
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
use vars qw/$opt_d $opt_o $opt_i $opt_r $opt_g $opt_c/;
getopts("d:r:g:i:o:c:");
my $usage =
".USAGE.   
rs_hisat.pl -d < paired end fastq.gz files > -i < file index number to process > -r < path to reference directory > -g < gtf/gff reference annotation > -o < output dir >

.DESCRIPTION.

.OPTIONS.
  -d  paired end fastq.gz files
  -i  index number of the file to process
  -r  path to reference transcriptome index (prefix of *.*ht2 files)
  -g  gtf/gff reference annotation
  -c  number of cores to use (optional - default 16)
  -o  output directory for alignment files
  
.KEYWORDS.
hisat, transcriptome, alignment, realignment, bam, fastq
\n";

die $usage unless defined $opt_d
              and defined $opt_i
              and defined $opt_r
              and defined $opt_g 
              and defined $opt_o;

my $FQ_DIR  = $opt_d;
my $INDEX   = $opt_i;
my $GTFPATH = $opt_g;
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
if (! -e "$OUT_DIR/hisat2"){
  `mkdir $OUT_DIR/hisat2`;
  `chmod ug+rw $OUT_DIR/hisat2`;
}
if (! -e "$OUT_DIR/stringtie"){
  `mkdir $OUT_DIR/stringtie`;
  `chmod ug+rw $OUT_DIR/stringtie`;
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
if (-e "$OUT_DIR/status/$id.hisat2.ck"){
  die "$id has already been processed according to its status file.\n";
}else{
  `echo processing > $OUT_DIR/status/$id.hisat2.ck`;
}

if (!defined($FQ_PAIRS{$id}{"R1"}) or !defined($FQ_PAIRS{$id}{"R2"})){
  die "Missing R1 or R2 pair for $id\n";
}
 
print "HISAT RUNS FOR $id ...\n"; 
`mkdir $OUT_DIR/hisat2/$id`;
`mkdir $OUT_DIR/stringtie/$id`;

print "HISAT RUNS FOR $id ...\n";
# HISAT2 alignment
print "hisat2 -x $REF_DIR -1 $FQ_PAIRS{$id}{R1} -2 $FQ_PAIRS{$id}{R2} -p $max_cores --dta -S $OUT_DIR/hisat2/$id/$id\_hisat2.sam";
      `hisat2 -x $REF_DIR -1 $FQ_PAIRS{$id}{R1} -2 $FQ_PAIRS{$id}{R2} -p $max_cores --dta -S $OUT_DIR/hisat2/$id/$id\_hisat2.sam`;

`/users/jrwhite/jlib/samtools-1.3/samtools view --threads $max_cores -Su $OUT_DIR/hisat2/$id/$id\_hisat2.sam  | /users/jrwhite/jlib/samtools-1.3/samtools sort --threads $max_cores -m 1000M -o $OUT_DIR/hisat2/$id/$id\_hisat2.sorted.bam -T $OUT_DIR/hisat2/$id/$id\_hisat2.tmp`;

print "stringtie...\n";
print "stringtie -B -G $GTFPATH -p $max_cores $OUT_DIR/hisat2/$id/$id\_hisat2.sorted.bam -o $OUT_DIR/stringtie/$id/$id\.gff -A $OUT_DIR/stringtie/$id/$id.gene_abund.txt\n";
      `stringtie -B -G $GTFPATH -p $max_cores $OUT_DIR/hisat2/$id/$id\_hisat2.sorted.bam -o $OUT_DIR/stringtie/$id/$id\.gff -A $OUT_DIR/stringtie/$id/$id.gene_abund.txt`;

