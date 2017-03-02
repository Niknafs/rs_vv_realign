#!/usr/bin/perl
##*********************************************************************
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_t $opt_n $opt_o $opt_c/;
getopts("t:n:o:c:");
my $usage =
".USAGE.   
rs_mito_TN_pipeline.pl -t < tumor fastq.gz file(s) > -n < normal fastq.gz file(s) > -o < output dir >

.DESCRIPTION.
Given a tumor and normal pair of fastq.gz files, this pipeline performs screening of human genome, 
followed by alignment to the human mitocondrial genome. 

.OPTIONS.
  -d  directory of paired end fastq.gz files
  -i  index number of the file to process
  -c  number of cores to use (optional - default 16)
  -o  output directory for results
  
.KEYWORDS.
\n";

die $usage unless defined $opt_d
              and defined $opt_i
              and defined $opt_o;

my $FQ_DIR  = $opt_d;
my $INDEX   = $opt_i;
my $OUT_DIR = $opt_o;
my $max_cores = 16;

my $PICARD_DIR    = "/users/jrwhite/jlib/picard-tools-1.84";
my $REF_DIR_MT    = "/dcl01/scharpf/data/reference/grch38/grch38/MT.fa";
my $REF_DIR_NONMT = "/dcl01/scharpf/data/reference/grch38/grch38/HomoSapiensGRCh38noMT.fa";

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
if (-e "$OUT_DIR/status/$id.bt2.ck"){
  die "$id has already been processed according to its status file.\n";
}else{
  `echo processing > $OUT_DIR/status/$id.bt2.ck`;
}

if (!defined($FQ_PAIRS{$id}{"R1"}) or !defined($FQ_PAIRS{$id}{"R2"})){
  die "Missing R1 or R2 pair for $id\n";
}
 
print "BT2 RUNS AGAINST HUMAN non MT FOR $id ...\n"; 
`mkdir $OUT_DIR/$id`;
`mkdir $OUT_DIR/$id/nonMT_screen`;
`mkdir $OUT_DIR/$id/MT_alignment`;

# now align to nonMT
`bowtie2 -x $REF_DIR_NONMT -1 $FQ_PAIRS{$id}{R1} -2 $FQ_PAIRS{$id}{R2} -p $max_cores --sensitive-local --no-unal --un-conc-gz $OUT_DIR/$id/nonMT_screen/$id\.bothunmapped.fastq.gz -S $OUT_DIR/$id/nonMT_screen/$id\.mapped.sam`;

# now align to the correct MT
`bowtie2 -x $REF_DIR_MT -1 $OUT_DIR/$id/nonMT_screen/$id\.bothunmapped.fastq.1.gz -2 $OUT_DIR/$id/nonMT_screen/$id\.bothunmapped.fastq.2.gz -p $max_cores --sensitive-local --no-unal --un-conc-gz $OUT_DIR/$id/MT_alignment/$id\.bothunmapped.fastq.gz -S $OUT_DIR/$id/MT_alignment/$id\.mapped.sam`;

# covert to sorted indexed bam file
`samtools view -bS $OUT_DIR/$id/MT_alignment/$id\.mapped.sam | samtools sort -O bam -o $OUT_DIR/$id/MT_alignment/cleaned.sorted.bam -T $OUT_DIR/$id/MT_alignment/cleaned`;
`samtools index $OUT_DIR/cleaned.sorted.bam`;

# mpileup and VarScan2 calls:

`samtools mpileup --fasta-ref $REF_DIR_MT $OUT_DIR/$id/MT_alignment/cleaned.sorted.bam > $OUT_DIR/$id/MT_alignment/cleaned.sorted.mpileup`;
# VARSCAN2 OPTIONS *********************************************************************
# --min-coverage	Minimum read depth at a position to make a call [8]
# --min-reads2	Minimum supporting reads at a position to call variants [2]
# --min-avg-qual	Minimum base quality at a position to count a read [15]
# --min-var-freq	Minimum variant allele frequency threshold [0.01]
# --min-freq-for-hom	Minimum frequency to call homozygote [0.75]
# --p-value	Default p-value threshold for calling variants [99e-02]
# --strand-filter	Ignore variants with >90% support on one strand [1]
# --output-vcf	If set to 1, outputs in VCF format
# --variants	Report only variant (SNP/indel) positions (mpileup2cns only) [0]
# 
my $vmincoverage   = 50;
my $vminreads2     = 4;
my $vminavgqual    = 25;
my $vminvarfreq    = 0.001;
my $vminfreqforhom = 0.75; 

`java -jar ~/jlib/VarScan.v2.4.2.jar mpileup2cns $OUT_DIR/$id/MT_alignment/cleaned.sorted.mpileup --min-coverage $vmincoverage --min-reads2 $vminreads2 --min-avg-qual $vminavgqual --min-var-freq $vminvarfreq --min-freq-for-hom $vminfreqforhom > $OUT_DIR/$id/MT_alignment/VarScan2.output.txt`;

# Mark Duplicates and remove
`java -Xms8g -Xmx32g -jar $PICARD_DIR/MarkDuplicates.jar I=$OUT_DIR/$id/MT_alignment/cleaned.sorted.bam O=$OUT_DIR/$id/MT_alignment/collapsed.bam METRICS_FILE=$OUT_DIR/$id/MT_alignment/collapsed.bam.metrics ASSUME_SORTED=True REMOVE_DUPLICATES=True CREATE_INDEX=True`;

# To extracted mapped reads where both mates mapped:
`samtools view -b -F12 $OUT_DIR/$id/MT_alignment/collapsed.bam > $OUT_DIR/$id/MT_alignment/collapsed.map-map.bam`;

# To extract mapped reads whose mates are unmapped:
`samtools view -b -F4 -f8 $OUT_DIR/$id/MT_alignment/collapsed.bam > $OUT_DIR/$id/MT_alignment/collapsed.map-unmap.bam`;

my $bed = "$OUT_DIR/$id/MT_alignment/MT.bed";
open OUT, ">$bed" or die;
print OUT "MT\t0\t16568\n";
close OUT;

# MPILEUP
open MM, "samtools mpileup -l $bed -d 1e6 -AB $OUT_DIR/$id/MT_alignment/collapsed.map-map.bam |";
open MMO, ">$OUT_DIR/$id/MT_alignment/map-map.coverage.txt";
print MMO "Chromosome\tCoord\tTotal\tA\tC\tG\tT\n";
while (<MM>){
  # chr10.fa        54860   N       2       gG      i~
  my @A = split "\t", $_;
  my $str = uc($A[4]);
  my @str = split "", $str;
  my %counts = (
    "A" => 0,
    "C" => 0,
    "G" => 0,
    "T" => 0
  );
  foreach my $s (@str){
    $counts{$s}++;
  }
  print MMO "$A[0]\t$A[1]\t$A[3]\t$counts{A}\t$counts{C}\t$counts{G}\t$counts{T}\n";
}
close MM;
close MMO;

open MU, "samtools mpileup -l $bed -d 1e6 -AB $OUT_DIR/$id/MT_alignment/collapsed.map-unmap.bam |";
open MUO, ">$OUT_DIR/$id/MT_alignment/map-unmap.coverage.txt";
print MUO "Chromosome\tCoord\tTotal\tA\tC\tG\tT\n";
while (<MU>){
  # chr10.fa        54860   N       2       gG      i~
  my @A = split "\t", $_;
  my $str = uc($A[4]);
  my @str = split "", $str;
  my %counts = (
    "A" => 0,
    "C" => 0,
    "G" => 0,
    "T" => 0
  );
  foreach my $s (@str){
    $counts{$s}++;
  }
  print MUO "$A[0]\t$A[1]\t$A[3]\t$counts{A}\t$counts{C}\t$counts{G}\t$counts{T}\n";
}
close MU;
close MUO;

# END


