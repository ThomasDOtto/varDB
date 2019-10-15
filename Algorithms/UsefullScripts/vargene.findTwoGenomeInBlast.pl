#! /usr/bin/perl -w
#
# File: vargene.findSequenceInBlast.pl
# Time-stamp: <15-Nov-2013 15:40:23 tdo>
# $Id: $
#
# Copyright (C) 2013 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#
use strict;


my $f=shift;


my $ref_l=getList($f);
my $ref_r=getList(shift);

while (<STDIN>) {
  my @ar=split(/\t/);
 my ($g1)= $ar[0] =~ /^(\S+)\.g/;
my ($g2)= $ar[1] =~ /^(\S+)\.g/;

  if (! defined($g1)){
	($g1)= $ar[0] =~ /^(\S+)_/;
	($g2)= $ar[1] =~ /^(\S+)_/;
}
  if (defined($$ref_l{$g1}) && defined($$ref_r{$g2})) {
	print 
  }
}

sub getList{
  open F, shift or die "gimme list\n";
  my %h;
  while (<F>) {
	chomp;
	$h{$_}=1
  }
  return \%h
}
