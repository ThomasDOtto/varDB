#! /usr/bin/perl -w
#
# File: vargene.getdomains.Matrix.pl
# Time-stamp: <02-Jun-2014 15:07:29 tdo>
# $Id: $
#
# Copyright (C) 2014 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#


my %h;
my %l;
my $lim=shift;

if (!defined($lim)) {
  $lim=1;
  
}

while (<STDIN>) {
  chomp;
  
  my @ar=split(/\t/);
  $ar[1] =~ s/-$//g;
  my @b=split(/-/,$ar[1]);
  
  for (my $i=0; $i<(scalar(@b)-1); $i++) {
	   $l{$b[$i]}++;
		$h{$b[$i]}{$b[($i+1)]}++;
		
  }
}



print "Domains";
foreach my $n (sort keys %l) {
  if ($l{$n} < $lim) {
#	delete ($l{$n})
  }else {
    print "\t$n";
  }
  
}
print "\n";

foreach my $n (sort keys %l) {
  print $n;
  
  foreach my $m (sort keys %l) {
   
	if (defined($h{$n}{$m})) {
	  print "\t$h{$n}{$m}";
	}
	else {
	  print "\t0"
	}
  }
print "\n";
}
