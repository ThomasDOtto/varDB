#! /usr/bin/perl -w
#
# File: vargene.hmm2gff.pl
# Time-stamp: <01-Oct-2015 11:59:27 tdo>
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
my $file=shift;

my $scorelimit=shift;
#my $limitfile = shift;

use Data::Dumper;

my ($ref_score, $ref_length) =getLimit($file);

if ((!defined($scorelimit))) {
  $scorelimit=10
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


### for overlap
my %overlap;

while (<STDIN>) {
  my @ar=split(/\s+/);
  ##                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
  # target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
  #------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
  #PA0032-C.g62         -           1413 ATS                  -            466         0 1132.5  95.1   1   2  1.9e-151  1.4e-150  503.5  22.2     4   293   764  1023   761  1033 0.95 -
  #  Pf3D7_01_v3     embl    gene    53169   53280   1000    -       .

  my $id=$ar[3];
  my $lgth=($ar[16]-$ar[15]);

  
  if (!defined($$ref_length{$ar[0]})) {
	print $ar[0]."\n"; print ; print $$ref_length{$ar[0]};
  }
  
 
  if ($ar[13]>=$scorelimit &&
	  $ar[13]>= (0.3 * $$ref_score{$ar[0]} ) ### 60% of maxscore
	  && $lgth > (0.8 *$$ref_length{$ar[0]} )
	 ) {
	my $col=1;
	my $forCol=$ar[0];
	$forCol =~ s/cluster\.\d+\.//g;
	
	my ($root) = $forCol;#$ar[0] =~ /(\S+)/;
	  
	if (defined($root) && ($color{$root})){
	  $col=$color{$root};
	} else {
	  $col=1
	}
	
	
	#### check here for overlap...
	
	
	my $ok=1;
	
	for (my $i=($ar[17]+20); $i<=($ar[18]-20);$i++) {
	  if (defined($overlap{$ar[3]}[$i])) {
		$ok=0;
		$i=($ar[18]+1)
	  }
	}

	if ($ok) {
	### now set the overlap to 1
	 for (my $i=($ar[17]+20); $i<=($ar[18]-20);$i++) { 
		$overlap{$ar[3]}[$i]=1;
	  }
	  
	}
	
	$ar[0] =~ s/Domain.//g;
	if ($ok) {
	#  print ">>$scorelimit $ar[13] $$ref_score{$ar[0]} $$ref_length{$ar[0]} $lgth\n";
	#  print "$ar[3]\tpfam\tCDS\t".(3*$ar[17]+1)."\t".(3*$ar[18]+1)."\t$ar[13]\t+\t.\tnote=$ar[0]-score=$ar[13];label=$ar[0];color=".$col."\n";
	  $h{$ar[3]}{$ar[17]}="$ar[3]\tpfam\tCDS\t".(3*$ar[17]+1)."\t".(3*$ar[18]+1)."\t$ar[13]\t+\t.\tnote=$ar[0]-score=$ar[13];label=$ar[0];color=".$col."\n";
	}
  }
}

foreach my $g (sort keys %h) {
  foreach my $p (sort {$a <=> $b}  keys %{$h{$g}} ) {
	print $h{$g}{$p};
  }
}



sub getLimit{

  open F, shift or die "Problems\n"; # "/lustre/scratch108/parasites/tdo/Pfalciparum/VAR/Assembly.Version2/Analysis/SevenVARgenes/Vardom.Domain.info.txt" or die "Problems...\n";

  my %length;
  my %score;
  while (<F>) {
	chomp;
	my ($id,$length,$score) = split(/\t/);
	if (!defined($score)) {
	  $score =10
	}
	$length{$id}=$length;
	$score{$id} = $score;
	
  }
  
  return (\%score,\%length)
}
