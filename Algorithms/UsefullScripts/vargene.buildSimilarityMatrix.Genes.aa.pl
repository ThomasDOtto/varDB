#! /usr/bin/perl -w
#
# File: vargene.coverageVarGenes.pl
# Time-stamp: <28-Jun-2016 10:54:32 tdo>
# $Id: $
#
# Copyright (C) 2012 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Descript
### sum = 2*alignemnts/( sum g1 + sum g2)
use strict;

my ($ref_fasta, $ref_length)=getFasta(shift);
my $result=shift;

my %seen;
my %h;
my %length;
my $MIN=0;



while (<STDIN>) {
# PA0050_VAR_12 PH0026_VAR_3 99.88 830 1 0 6704 7533 4414 5243 0.0
# 1637 7751 5414
  chomp;
  
  my ($q,$s,$id,$ov,$one,$two,$qStart,$qEnd,$sStart,$sEnd,
	  $score,$value,$qLength,$sLength)=split(/\t/);
  
 if ($sEnd<$sStart) {
	my $tmp=$sEnd;$sEnd=$sStart;$sStart=$tmp;	
  }

  
  my ($sGenome) = $s;# =~ /^(\S+)\.g/;
  my ($qGenome) = $q;# =~ /^(\S+)\.g/;

  
  
  if (!defined($sGenome)){
	
	($sGenome) = $s;# =~ /^(\w+)-/;
  } 
  if (!defined($qGenome)){
	($qGenome) = $q;# =~ /^(\w+)-/;
  }
  if (!defined($qGenome)){ print $q };
  if (!defined($sGenome)){ print $s };
  
#  if (!defined($seen{$s})) {
#  	$seen{$s}=1;
#  }
  
  
  ### carefull notSeen will read %seen
  if (!(defined($seen{$q}{$sGenome}))){#} && $sGenome ne $qGenome) {
	### befor in the if  && notSeen($q,$sGenome,$qStart,$qEnd
	$seen{$q}{$sGenome}=1;
	$h{$qGenome}{$sGenome}+=($id*$ov/100);
	#	for my $i ($qStart..$qEnd){
	#  if (($i%20)==0){
	#$seen{$q}{$sGenome}{$i}=1;	
	#  }
	#}
  }
  
  
}

my %normed;
foreach my $query (keys %$ref_length) {
  foreach my $subject (keys %$ref_length){
	if (defined($h{$query}{$subject})){
	  if (!defined($h{$subject}{$query})) {
		$h{$subject}{$query}=0
	  }
	  $normed{$query}{$subject}=(100*( ($h{$query}{$subject}+$h{$subject}{$query})/($$ref_length{$subject}+$$ref_length{$query})));
	  if ($normed{$query}{$subject}>100) {
		$normed{$query}{$subject}=100
	  }
	  if ($normed{$query}{$subject}<$MIN) {
		$normed{$query}{$subject}=$MIN
	  }
	}
	else {
	  $normed{$query}{$subject}=$MIN;
	}
  }  
}

open F, "> $result" or die "Probelm\n";
print F "Genome";
foreach my $query (sort keys %$ref_length) {
	print F "\t$query";
}
print F "\n";
foreach my $query (sort keys %$ref_length) {
	print F $query;
	foreach my $subject (sort keys %$ref_length){
	   print F "\t$normed{$query}{$subject}";
	}
	print F "\n";
}
close(F);

sub notSeen{
	my $q= shift;
	my $genome=shift;
	my $s=shift;
	my $e=shift;
	
	for ($s..$e){
		if (defined($seen{$q}{$genome}{$_})){
			return 0	
		}	
	}	
	
	return 1;	
}

sub min{
  my $v1=shift;
  my $v2=shift;
  if ($v1< $v2) {
	return $v1
  }
  else {
	return $v2
  }
  
}
sub getFasta{
  my $fasta=shift;
  
  open F, $fasta or die "prob couldn't find fasta file: $fasta  \n";
  
  my %h;
  my %length;
  
  my $name;
  
  while (<F>) {
        chomp;
        if (/^>(\S+)/){
		  my $s=$1;
		  
		  ($name) = $s;# =~ /^(\S+)\.g/;
        }
        else {
          chomp;
		  $length{$name} += length($_);
		  
          $h{$name}.=$_;
        }
  }
  return (\%h,\%length);
}
