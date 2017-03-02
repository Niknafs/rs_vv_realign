#!/usr/bin/perl
##*********************************************************************
##  rs_project_setup.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##  created: 2015-06-18
##*********************************************************************
## This script creates the standardized directory structure necessary for a series
# of human genome analyses
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
rs_project_setup.pl -d < name of project directory >

.OPTIONS.
  -d  project directory name (e.g. CRC_EXOME)
\n";

die $usage unless defined $opt_d;

my $DIR  = $opt_d;

my @aligners = qw/eland novoalign/;
my @refs     = qw/hg18 hg19 GRCh38/; 

`mkdir $DIR`;
`mkdir $DIR/metadata`;
`touch $DIR/DESCRIPTION.TXT`;
`mkdir $DIR/fastq`;
`mkdir $DIR/origin_data`;
`touch $DIR/origin_data/README`;

foreach my $r (@refs){
  `mkdir $DIR/$r`;
foreach my $a (@aligners){
  `mkdir $DIR/$r/$a`;
  `mkdir $DIR/$r/$a/bam`;
  `mkdir $DIR/$r/$a/doc`;
}
}

