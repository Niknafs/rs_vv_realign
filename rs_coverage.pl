#!/usr/bin/perl
#*********************************************************************
#  rs_coverage.pl*
#  Author:  James Robert White, PhD
#  Email:   james.dna.white@gmail.com
#  Created: 2015-10-01
#*********************************************************************

#*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use File::Spec;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/;
#*********************************************************************

use vars qw/$opt_o $opt_i/;
getopts("i:o:");
my $usage =
".USAGE.
rs_coverage.pl -i < sorted indexed bam file > -o < output directory >

TO RUN INTERACTIVELY YOU MUST REQUEST A LARGE MEM_FREE RESOURCE E.G.
qrsh -l mf=40G,h_vmem=40G -l gwas -q gwas.q\@compute-059 -pe local 1

NOTE THE INPUT BAM FILE MUST BE SORTED AND INDEXED

.KEYWORDS.
coverage calculation
\n";

die $usage unless defined $opt_i
              and defined $opt_o;

my $bam     = File::Spec->rel2abs( $opt_i );
my $outDir  = $opt_o;

print "Input bam: $bam\n";
print "Specified output directory: $outDir\n";

# PROCESSING ****************************************************************
# PARAMETERS
my $PICARD_DIR           = "/users/jrwhite/jlib/picard-tools-1.84";

if (! -e "$outDir"){
  `mkdir $outDir`;
  `chmod ug+rw $outDir`;
}

# Begin by reheadering the bam artificially increasing the header lines
`samtools view -H $bam > $outDir/header.sam`;
open IN, "$outDir/header.sam" or die;
open OUT, ">$outDir/reheader.sam" or die;
while(<IN>){
  chomp($_);
  if ($_ =~ /^\@SQ/){
    my @A = split /\tLN\:/, $_;
    $A[1] = $A[1]+1000;
    $_ = join("\tLN:", @A);
  }
  print OUT "$_\n";
}
close IN;
close OUT;

# REHEADER
`samtools reheader $outDir/reheader.sam $bam > $outDir/orig_reheader.bam`;

# Clean up a messy bam file
`java -Xms8g -Xmx32g -XX:ParallelGCThreads=1 -jar $PICARD_DIR/CleanSam.jar I=$outDir/orig_reheader.bam O=$outDir/cleaned.bam`;

# Mark Duplicates and remove
`java -Xms8g -Xmx32g -jar $PICARD_DIR/MarkDuplicates.jar I=$outDir/cleaned.bam O=$outDir/collapsed.bam METRICS_FILE=$outDir/collapsed.bam.metrics ASSUME_SORTED=True REMOVE_DUPLICATES=True CREATE_INDEX=True`;

# To extracted mapped reads where both mates mapped:
`samtools view -b -F12 $outDir/collapsed.bam > $outDir/collapsed.map-map.bam`;

# To extract mapped reads whose mates are unmapped:
`samtools view -b -F4 -f8 $outDir/collapsed.bam > $outDir/collapsed.map-unmap.bam`;

# Create counts:
`bedtools genomecov -bga -ibam $outDir/collapsed.map-unmap.bam > $outDir/map-unmap.coverage.txt`;
`bedtools genomecov -bga -ibam $outDir/collapsed.map-map.bam   > $outDir/map-map.coverage.txt`;

# END
