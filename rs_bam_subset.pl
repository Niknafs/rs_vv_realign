#!/usr/bin/perl
#*********************************************************************
#  rs_bam_subset.pl*
#  Author:  James Robert White, PhD
#  Email:   james.dna.white@gmail.com
#  Created: 2016-02-04
#*********************************************************************

#*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use File::Spec;
use warnings;
use strict;
use Switch;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/;
#*********************************************************************

use vars qw/$opt_o $opt_b $opt_i/;
getopts("i:b:o:");
my $usage =
".USAGE.
rs_bam_subset.pl -i < input bam file > -b < input bed file of regions > -o < output subsetted bam file >

.PARAMETERS.
-i   input bam file 
-b   input bed file of regions (3-column format e.g. chrX<tab>000000<tab>000000<newline>
-o   output subsetted bam file -- must end in .bam

TO RUN INTERACTIVELY YOU MUST REQUEST A LARGE MEM_FREE RESOURCE E.G.
qrsh -l mf=40G,h_vmem=40G -l gwas -q gwas.q\@compute-059 -pe local 1

NOTE THE FOLLOWING MODULES ARE REQUIRED
module load perl
module load samtools/1.1
module load bedtools

.KEYWORDS.
coverage calculation
\n";

die $usage unless defined $opt_i
              and defined $opt_b
              and defined $opt_o;

my $bam        = File::Spec->rel2abs( $opt_i );
my $bed        = $opt_b;
my $outbam     = $opt_o;
my $PICARD_DIR = "/users/jrwhite/jlib/picard-tools-1.84";

if ($outbam !~ /.bam$/){
  die "ERROR: Output bam file doesnt end in .bam\n\n";
}

if (-e $outbam){
  die "ERROR: the output bam file $outbam already exists!\n";
}

if (! -e $bam){
  die "ERROR: the input bam file $bam does not exist!\n";
}

if (! -e $bed){
  die "ERROR: the input bam file $bed does not exist!\n";
}

print "*********************************************************\n";
print "Input bam file: $bam\n";
print "Input bed file: $bed\n";
print "Output bam: $outbam\n";
print "*********************************************************\n";

# Check bam file
print "Checking bam file header for seq identifiers...\n";
# @SQ	SN:chrM	LN:16571
#samtools view -H LP6008081-DNA_B02.merged.bam | grep "^@SQ"
my $sqs = `samtools view -H $bam | grep \"^\@SQ\"`;
chomp($sqs);
my @sqs = split "\n", $sqs;

my %seqids = ();
foreach my $sq (@sqs){
  my @sq = split "\t", $sq;
  print "$sq\n";
  $sq[1] =~ s/SN://g; 
  $seqids{$sq[1]} = 1;
}

# Check BED file
print "Checking bed file for consistency...\n";
open IN, "$bed" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if ($#A > 2){
    die "Extra spaces found in the bed file $bed\n";
  }

  if (!defined($seqids{$A[0]})){
    die "ERROR: Seq identifier conflict. $A[0] is not found in the bam file header!\n\n";
  } 
}
close IN;

print "Bam and Bed files appear consistent in sequence identifiers.\n";
print "*********************************************************\n";
# BEGIN REGION FILTERING

print "Identifying reads that hit the ROIs...\n";  
`samtools view -b -L $bed $bam > $outbam.mapped.bam`;
`samtools view $outbam.mapped.bam | cut -f 1 | sort -u > $outbam.mapped.txt`;

print "Collecting unmapped and offmapped read pairs...\n";
`java -Xms8g -Xmx32g -jar $PICARD_DIR/FilterSamReads.jar I=$bam O=$outbam FILTER=includeReadList READ_LIST_FILE=$outbam.mapped.txt WRITE_READS_FILES=FALSE`;

`rm $outbam.mapped.txt $outbam.mapped.bam`;

# END
