#! /usr/bin/perl -w
#
# File: vargene.parseOrthomcl.pl
# Time-stamp: <06-Feb-2013 17:16:04 tdo>
# $Id: $
#
# Copyright (C) 2013 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#

my %con;

#$con{PP}="Peru";
$con{PF}="Ghana";
$con{PA}="Gambia";
$con{PC}="Kenia";
$con{PM}="Mali";
$con{PT}="Malauwi";
#$con{PK}="Burkina Faso";
#$con{PE}="Tansamina";
$con{PS}="Senegal";
#$con{PR}="Bangladesh";
$con{PD}="Thailand";
$con{PH}="Cambodia";
#$con{PN}="PNG";
$con{PV}="Vietnam";
#$con{PZ}="PZ";
#$con{PX}="Control";
$con{QG}="Congo";
$con{QE}="Laos";
$con{PU}="Guinea";
#$con{PL}="Different 3D7 and mixed";
#$con{PX}="Sanger Control";

my @order=("PS","PA","PU","PM","PF","QG","PT","PC","PD","QE","PH","PV");
print "ortho\tgenes\tgenomes\t";

#  foreach my $k (keys %con) {
foreach my $k (@order){
	print "$con{$k}\t";
  }
print "\n";

while (<STDIN>) {
  chomp;

  my @ar=split(/\s+/);
  # ORTHOMCL8(50 genes,50 taxa):
  
  my ($n,$g)=$ar[0]=~/(ORTHOMCL\d+)\((\d+)/;
  
  my ($t) = $ar[1] =~ /genes,(\d+)/;
  
  $ar[0]='(XX';
  my %h;
  foreach (@ar) {
	/\((\S{2})/;
	$h{$1}++;
  }
  print "$n\t$g\t$t\t";

foreach my $k (@order){  
#  foreach my $k (keys %con) {
	if (defined($h{$k})) {
	  print "$h{$k}\t";
	}
	else {
	  print "0\t"
	}
  }
  print "\n";
  
}
