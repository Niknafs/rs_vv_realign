#!/usr/bin/perl
##*********************************************************************
##  rs_xenotrans.pl
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2016-06-02
##*********************************************************************
## This script performs processing of fastq.gz files
## 
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
rs_xenotrans.pl -d < directory of paired end fastq.gz files > -i < file index number to process > -r < path to human reference directory > -o < output dir >

.DESCRIPTION.
This pipeline performs automated transcriptome quantification using STAR/RSEM for 
xenograft samples. Mouse (v. GRCm38) is also searched for and cross-species matches are removed
from consideration. 

.OPTIONS.
  -d  paired end fastq.gz files
  -i  index number of the file to process
  -r  path to reference transcriptome index (prefix of RSEM-STAR index files)
      e.g. /dcl01/scharpf/data/reference/grch37/starrsem_GRCh37_noMTrRNA 
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

my $MOUSE_REF_DIR = "/dcl01/scharpf/data/reference/Mmusculus_mm10/grcm38_snp_tran/genome_snp_tran";

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

# status rsem 
 
`mkdir $OUT_DIR/rsem/$id`;
`mkdir $OUT_DIR/rsem/$id/contaminant_ck`;
print "HISAT RUNS FOR $id ...\n";
# HISAT2 alignment
`hisat2 -x $MOUSE_REF_DIR -1 $FQ_PAIRS{$id}{R1} -2 $FQ_PAIRS{$id}{R2} -p $max_cores --no-unal --un-conc-gz $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.fastq.gz --dta -S $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.mapped.sam`;

`mv $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.fastq.1.gz $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.1.fastq.gz`;
`mv $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.fastq.2.gz $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.2.fastq.gz`;

# $id_hisat2.bothunmapped.fastq.1.gz
 print "STAR-RSEM RUNS AGAINST HUMAN (NO Mitochondria or rRNA) FOR $id ...\n";
 `mkdir $OUT_DIR/rsem/$id/filtered`;
 `rsem-calculate-expression --paired-end --star --star-path /users/jrwhite/jbin --star-gzipped-read-file -p $max_cores $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.1.fastq.gz $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.2.fastq.gz $REF_DIR $OUT_DIR/rsem/$id/filtered/$id`;

# fusions!
  print "STAR-RSEM FUSION DETECTION for $id...\n";
  `/users/jrwhite/jlib/STAR-Fusion/STAR-Fusion --left_fq $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.1.fastq.gz --right_fq  $OUT_DIR/rsem/$id/contaminant_ck/$id\_hisat2.bothunmapped.2.fastq.gz -O $OUT_DIR/rsem/$id/filtered/star_fusion_outdir --genome_lib_dir /dcl01/scharpf/data/reference/grch37/GRCh37_gencode_v19_CTAT_lib --verbose_level 2`;

print "COMPLETED.\n";
