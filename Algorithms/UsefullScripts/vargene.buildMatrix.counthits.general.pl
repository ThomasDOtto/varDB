#! /usr/bin/perl -w
#
# File: vargene.coverageVarGenes.pl
# Time-stamp: <03-Jun-2014 10:09:08 tdo>
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


while (<STDIN>) {
# PA0050_VAR_12 PH0026_VAR_3 99.88 830 1 0 6704 7533 4414 5243 0.0
# 1637 7751 5414
  chomp;
  
  my ($q,$s,$id,$ov,$one,$two,$qStart,$qEnd,$sStart,$sEnd,
	  $score,$value,$qLength,$sLength)=split(/\t/);
  
 if ($sEnd<$sStart) {
	my $tmp=$sEnd;$sEnd=$sStart;$sStart=$tmp;	
  }

  
  my ($sGenome) = $s =~ /^(\S+)_/;
  my ($qGenome) = $q =~ /^(\S+)_/;

  
  
  if (!defined($sGenome)){
	
	($sGenome) = $s =~ /^(\w+)-/;
  } 
  if (!defined($qGenome)){
	($qGenome) = $q =~ /^(\w+)-/;
  }
  if (!defined($qGenome)){ print $q };
  if (!defined($sGenome)){ print $s };
  
  if (!defined($seen{$q}{$sGenome})) {
 	$h{$qGenome}{$sGenome}++;
  	$seen{$q}{$sGenome}=1;
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
	  if (defined($h{$query}{$subject})) {
		print F "\t$h{$query}{$subject}";
	  }
	  else {
		print F "\t0"
	  }
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
		  
		  ($name) = $s =~ /^(\S+)_/;
        }
        else {
          chomp;
		  $length{$name} += length($_);
		  
          $h{$name}.=$_;
        }
  }
  return (\%h,\%length);
}
