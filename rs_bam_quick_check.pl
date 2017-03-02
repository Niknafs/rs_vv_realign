#!/usr/bin/perl
##*********************************************************************
##  rs_bam_to_fq.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-04-02
##*********************************************************************
## This script performs processing of bam alignments files including
#  .sorting 
#  .converting to fastq and formatting
#  .splitting for size
#   
##*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_d/;
getopts("d:");
my $usage =
".USAGE.   
rs_bam_quick_check.pl -d < directory containing bam files > 

.DESCRIPTION.
  Reports on the quality of all detected bam files in a directory 
  (and subdirectories). Describes files that it finds an issue with.
  Specifically looks for truncated files.

.OPTIONS.
  -d  dir of bam files
      ( can be complex with subdirectories 
        containing bam files as well 
        e.g. dir1/subdir2/subsubdir3/sample.bam )

.KEYWORDS.
\n";

die $usage unless defined $opt_d;

my $SAMTOOLS_PATH = "/users/jrwhite/jlib/samtools-1.3/samtools";
my $BAM_DIR = $opt_d;

if (! -e "$BAM_DIR"){
  die "Unable to locate input directory: $BAM_DIR\n";
}


my %BAM_FILES = ();
my $BAM_LIST  = `find $BAM_DIR -name "*.bam" -print`;
chomp($BAM_LIST);
my @BAM_LIST   = split "\n", $BAM_LIST;
my $total_bams = 0;
foreach my $bam (@BAM_LIST){
  $BAM_FILES{$bam} = 1;
  $total_bams++;
}

print "QUICK CHECK REPORT:\n";
print "Total bam files: $total_bams\n";
if ($total_bams == 0){
  die "No bams found.\n";
}
print "Beginning samtools (v1.3) quick checks...\n";
my $badcount = 0;
my $iter = 1;
foreach my $b (sort keys %BAM_FILES){
  my $prctdone = sprintf("%2.1f", 100*$iter/$total_bams);
  print "$iter ($prctdone\%)\r";
  $iter++;
  
  my $output = `$SAMTOOLS_PATH quickcheck -vvv $b 2>&1`;
  my $size   = `ls -lht $b | awk \'\{print \$5\}\'`;
  chomp($size);

  if ($output !~ /good EOF/){
    print "It appears $b is truncated...filesize=$size\n";
    $badcount++;
  }
}

print "\n*********\nSummary: $badcount bam file(s) showed some problems and may require rerunning.\n";
