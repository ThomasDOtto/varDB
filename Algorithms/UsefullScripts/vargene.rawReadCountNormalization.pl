#! /usr/bin/perl -w
#
# File: vargene.rawReadCountNormalization.pl
# Time-stamp: <23-Jul-2015 11:31:14 tdo>
# $Id: $
#
# Copyright (C) 2015 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#

use strict;


my $NORM=10000;
my $file_length=shift;
my $factor=shift;

my %h;

### get read length
open F, $file_length or die "problem $file_length\n";
while (<F>) {
  my @ar=split(/\t/);
  $h{$ar[0]}=$ar[1];
}
close(F);


while (<STDIN>) {
  chomp;
  my ($name,$c) = split(/\t/);
  ### thousands are for the gene length...
  my $res=($c*$NORM*1000/($h{$name}*$factor));
  printf ("$name\t%0.2f\n",$res);
  
}

