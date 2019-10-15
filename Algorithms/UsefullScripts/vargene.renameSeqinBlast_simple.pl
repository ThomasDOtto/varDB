#! /usr/bin/perl -w
#
# File: vargene.findSequenceInBlast.pl
# Time-stamp: <25-Feb-2013 17:19:57 tdo>
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


while (<STDIN>) {
  my @ar=split(/\t/);

  if (defined($$ref_l{$ar[0]})){
	        $ar[0]=$$ref_l{$ar[0]};
} 
if ( defined($$ref_l{$ar[1]})) {
#	$ar[0]=$$ref_l{$ar[0]};
	$ar[1]=$$ref_l{$ar[1]};
}
print join ("\t",@ar);
  
}

sub getList{
  open F, shift or die "gimme list\n";
  my %h;
 my %s;
  while (<F>) {
	chomp;
	my @ar=split(/\t/);
	$h{$ar[0]}=$ar[1];
  }
  return \%h
}
