#! /usr/bin/perl -w
#
# File: vargene.hmm2gff.pl
# Time-stamp: <15-Nov-2013 14:38:24 tdo>
# $Id: $
#
# Copyright (C) 2013 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#

use strict;


my %h;
my %color;
my $scorelimit=shift;
#my $limitfile = shift;

use Data::Dumper;

my ($ref_score, $ref_length) =getLimit();


if ((!defined($scorelimit))) {
  $scorelimit=200
}
$color{"ATS"}=1;
$color{"CIDRa"}=2;
$color{"CIDRb"}=3;
$color{"CIDRd"}=4;
$color{"CIDRg"}=5;
$color{"CIDRpam"}=6;
$color{"DBLa"}=7;
$color{"DBLb"}=8;
$color{"DBLd"}=9;
$color{"DBLe"}=10;
$color{"DBLg"}=11;
$color{"DBLpam1"}=12;
$color{"DBLpam2"}=13;
$color{"DBLpam3"}=14;
$color{"DBLz"}=15;
$color{"Duffy_binding_like"}=1;

while (<STDIN>) {
  my @ar=split(/\s+/);
  ##                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
  # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
  #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
  #PA0032-C.g62         -           1413 ATS                  -            466         0 1132.5  95.1   1   2  1.9e-151  1.4e-150  503.5  22.2     4   293   764  1023   761  1033 0.95 -
  #  Pf3D7_01_v3     embl    gene    53169   53280   1000    -       .

  my $id=$ar[0];
  my $lgth=($ar[16]-$ar[15]);
  
  if ($ar[13]>=$scorelimit &&
	  $ar[13]>= (0.3 * $$ref_score{$ar[3]} ) ### 60% of maxscore
	  && $lgth > (0.9 *$$ref_length{$ar[3]} )
	 ) {
	my $col=1;
	if (defined($color{$ar[3]})){
			  $col=$color{$ar[3]};
	
	}
	$h{$ar[0]}{$ar[17]}="$ar[0]\tpfam\tCDS\t".(3*$ar[17]+1)."\t".(3*$ar[18]+1)."\t$ar[13]\t+\t.\tnote=\"$ar[3]-score=$ar[13]\";label=\"$ar[3]\";color=".$col."\n";
  }
}

foreach my $g (sort keys %h) {
  foreach my $p (sort {$a <=> $b}  keys %{$h{$g}} ) {
	print $h{$g}{$p};
  }
}



sub getLimit{

#  open F, "/lustre/scratch108/parasites/tdo/Pfalciparum/VAR/Assembly.Version2/Analysis/SevenVARgenes/Vardom.Domain.info.txt" or die "Problems...\n";
open F, "/nfs/pathogen003/tdo/Pfalciparum/VAR/Assembly.Version2/Analysis/SevenVARgenes/Vardom.Domain.info.txt" or die "problems \n";
  my %length;
  my %score;
  while (<F>) {
	chomp;
	my ($id,$length,$score) = split(/\t/);
	$length{$id}=$length;
	$score{$id} = $score;
	
  }
  
  return (\%score,\%length)
}
