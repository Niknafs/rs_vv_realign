#!/usr/bin/perl
##*********************************************************************
##  rs_vvpodsync.pl*
##  author: james robert white, phd. 
##  email: james.dna.white@gmail.com
##*********************************************************************
## This script performs the creation of two scripts for SGE rsyncs between
## JHPCE and the POD in Victor Velculescu's lab
###*********************************************************************
use Data::Dumper;
use Getopt::Std;
use Getopt::Long;
use warnings;
use strict;
use POSIX qw/ceil floor/;
use List::Util qw/first max maxstr min minstr reduce shuffle sum/; 
##*********************************************************************
use vars qw/$opt_f $opt_t $opt_p/;
getopts("f:t:p:");
my $usage =
".USAGE.   
rs_vvpodsync.pl -f < from this pod directory > -t < to this directory > -p < prefix >

e.g.
rs_vvpodsync.pl -f '/data1/Data\\ from\\ PGDx/Breast/Feb\\ 2015\\ Samples/ELAND\\ Bam\\ files/PGDX395*WGS.bam*' -t 'BRCA_FEB_bam' -p brca_rsync


.DESCRIPTION.
  PERFORM AN RSYNC WITH THE POD

.OPTIONS.
  -f  *from* this pod directory (in single quotes with \\'s before spaces!)
      can be a path e.g. 
      '/data1/Data\\ from\\ PGDx/Breast/Feb\\ 2015\\ Samples/ELAND\\ Bam\\ files/'  
      or a path with regular expression e.g.
      '/data1/Data\\ from\\ PGDx/Breast/Feb\\ 2015\\ Samples/ELAND\\ Bam\\ files/PGDX395*WGS.bam*'

  -t  *to* this directory (same formatting style as -f above) e.g.
      'myproject/bam'

  -p  output script prefix

      <prefix>.dryrun.sh -- a dry run script you can execute to confirm 
                            the dirs are correctly set up
                            execute on the command-line as: 
                            ./<prefix>.dryrun.sh

      <prefix>.sge.sh    -- the real script to submit to SGE via qsub e.g.
                            qsub <prefix>.sge.sh 
       
      Note: if the dry-run fails on the command line please recheck and confirm
      the path of the *from* directory 
\n";

die $usage unless defined $opt_f 
              and defined $opt_t
              and defined $opt_p;

my $prefix = $opt_p;
my $from   = $opt_f; 
my $to     = $opt_t;

open OUT, ">$prefix.dryrun.sh" or die "\nError: Can't open $prefix.dryrun.sh!\n";
my $line = "rsync --dry-run -av vvpoduser\@10.99.4.253:\'$from\' $to";
print OUT "$line\n";
close OUT;

open OUT, ">$prefix.sge.sh" or die "\nError: Can't open $prefix.sge.sh!\n";
print OUT 
"\#\!/bin/bash
\#\$ -cwd
\#\$ -j y
\#\$ -l mem_free=4G
\#\$ -l h_vmem=24G
\#\$ -l h_fsize=100G
\#\$ -l h_rt=48\:00\:00
\#\$ -l rnet
rsync -av vvpoduser\@10.99.4.253:'$from' $to
";
close OUT;

`chmod a+x $prefix*.sh`;

