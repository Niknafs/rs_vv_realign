#!/usr/bin/perl
#*********************************************************************
#  rs_mutect_pipeline.pl*
#  Author:  James Robert White, PhD
#  Email:   james.dna.white@gmail.com
#  Created: 2015-04-09
#*********************************************************************
# PROTOTYPE QUALITY

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

use vars qw/$opt_t $opt_n $opt_o $opt_i/;
getopts("t:n:o:i:");
my $usage =
".USAGE.
 rs_mutect_pipeline.pl -t <tumor bam> -n <normal bam> -o <output directory> -i <interval list>

.DESCRIPTION.

.OPTIONS.

.KEYWORDS.
\n";

die $usage unless defined $opt_t
              and defined $opt_n
              and defined $opt_o
              and defined $opt_i;

my $tbam     = File::Spec->rel2abs( $opt_t );
my $nbam     = File::Spec->rel2abs( $opt_n );
my $outDir   = $opt_o;
my $INTERVAL_LIST = File::Spec->rel2abs( $opt_i );

print "Input tumor bam: $tbam\n";
print "Input normal ba: $nbam\n";
print "Specified output directory: $outDir\n";

# PROCESSING ****************************************************************
# PARAMETERS
my $LOW_MEM_NPROC        = 8;
my $HIGH_MEM_NPROC       = 2;
my $DB_DIR               = "/dcl01/scharpf/data/reference/with_alternates/hg19_wUns";
my $HUMAN_GENOME_REF     = "$DB_DIR/hg19_wfas.fa";
my $COSMIC_REF           = "b37_cosmic_v54_120711.mod.vcf";
my $DBSNP_REF            = "dbsnp_138.b37.vcf";

my $PICARD_DIR           = "/users/jrwhite/jlib/picard-tools-1.84";
my $MUTECT_DIR           = "/users/jrwhite/jlib/mutect-1.1.4";
my $JAVA6                = "/usr/lib/jvm/java-1.6.0/bin/java";

`mkdir $outDir`;

# Step 2.1
if (! -e "$outDir/tumor.filtered.bam"){
  `samtools view -bF 4 -q 0 $tbam > $outDir/tumor.filtered.bam`;
  `samtools view -bF 4 -q 0 $nbam > $outDir/normal.filtered.bam`;
}

# Step 2.2
if (! -e "$outDir/tumor.filtered.sortedfixed.bam"){
  `java -Xmx20g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=$outDir/tmpA25 -jar $PICARD_DIR/FixMateInformation.jar I=$outDir/tumor.filtered.bam O=$outDir/tumor.filtered.sortedfixed.bam MAX_RECORDS_IN_RAM=2000000 TMP_DIR=$outDir/tmpA25 SO=coordinate VALIDATION_STRINGENCY=STRICT`;
  `java -Xmx20g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=$outDir/tmpA25 -jar $PICARD_DIR/FixMateInformation.jar I=$outDir/normal.filtered.bam O=$outDir/normal.filtered.sortedfixed.bam MAX_RECORDS_IN_RAM=2000000 TMP_DIR=$outDir/tmpA25 SO=coordinate VALIDATION_STRINGENCY=STRICT`;
}

# Step 2.3
if (! -e "$outDir/tumor.filtered.sortedfixed.reorder.bam"){
`java -Xmx24g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=$outDir/tmpA25A -jar $PICARD_DIR/ReorderSam.jar I=$outDir/tumor.filtered.sortedfixed.bam O=$outDir/tumor.filtered.sortedfixed.reorder.bam R=$HUMAN_GENOME_REF TMP_DIR=$outDir/tmpA25A VALIDATION_STRINGENCY=STRICT`;
`java -Xmx24g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=$outDir/tmpA25A -jar $PICARD_DIR/ReorderSam.jar I=$outDir/normal.filtered.sortedfixed.bam O=$outDir/normal.filtered.sortedfixed.reorder.bam R=$HUMAN_GENOME_REF TMP_DIR=$outDir/tmpA25A VALIDATION_STRINGENCY=STRICT`;
}

# Step 2.4
if (! -e "$outDir/tumor.clean.bam"){
`java -Xmx24g -XX:ParallelGCThreads=8 -jar $PICARD_DIR/CleanSam.jar I=$outDir/tumor.filtered.sortedfixed.reorder.bam O=$outDir/tumor.clean.bam`;
`java -Xmx24g -XX:ParallelGCThreads=8 -jar $PICARD_DIR/CleanSam.jar I=$outDir/normal.filtered.sortedfixed.reorder.bam O=$outDir/normal.clean.bam`;
}

# Step 2.5
if (! -e "$outDir/tumor.sorted.collapsed.bam"){
`java -jar $PICARD_DIR/MarkDuplicates.jar I=$outDir/tumor.clean.bam O=$outDir/tumor.sorted.collapsed.bam METRICS_FILE=$outDir/tumor.sorted.collapsed.bam.metrics ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true`;
`java -jar $PICARD_DIR/MarkDuplicates.jar I=$outDir/normal.clean.bam O=$outDir/normal.sorted.collapsed.bam METRICS_FILE=$outDir/normal.sorted.collapsed.bam.metrics ASSUME_SORTED=true REMOVE_DUPLICATES=true CREATE_INDEX=true`;
}

# Step 2.6
if (! -e "$outDir/tumor.sorted.collapsed_addrg.bam"){
`java -Xmx24g -XX:ParallelGCThreads=8 -jar $PICARD_DIR/AddOrReplaceReadGroups.jar I=$outDir/tumor.sorted.collapsed.bam O=$outDir/tumor.sorted.collapsed_addrg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1`;
`java -Xmx24g -XX:ParallelGCThreads=8 -jar $PICARD_DIR/AddOrReplaceReadGroups.jar I=$outDir/normal.sorted.collapsed.bam O=$outDir/normal.sorted.collapsed_addrg.bam RGID=group2 RGLB=lib2 RGPL=illumina RGPU=unit2 RGSM=sample2`;

# Step 2.8
`java -Xmx24g -XX:ParallelGCThreads=8  -jar $PICARD_DIR/BuildBamIndex.jar I=$outDir/tumor.sorted.collapsed_addrg.bam`;
`java -Xmx24g -XX:ParallelGCThreads=8  -jar $PICARD_DIR/BuildBamIndex.jar I=$outDir/normal.sorted.collapsed_addrg.bam`;

}

# Step 2.10
if (! -e "$outDir/mutect.vcf"){
`$JAVA6 -Xmx24g -Djava.io.tmpdir=$outDir/tmp -jar $MUTECT_DIR/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $HUMAN_GENOME_REF --input_file:tumor $outDir/tumor.sorted.collapsed_addrg.bam --input_file:normal $outDir/normal.sorted.collapsed_addrg.bam --out $outDir/mutect.call_stats.txt --coverage_file $outDir/mutect.coverage.wig.txt --vcf $outDir/mutect.vcf --num_threads 8 --downsample_to_coverage 100000 --initial_tumor_lod 3.0 --tumor_lod 4.0 --fraction_contamination 0 --required_maximum_alt_allele_mapping_quality_score -10 --intervals $INTERVAL_LIST --enable_extended_output`;
}
