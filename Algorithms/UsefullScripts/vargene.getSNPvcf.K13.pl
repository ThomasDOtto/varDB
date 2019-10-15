#! /usr/bin/perl -w
#
# File: vargene.getSNPvcf.pl
# Time-stamp: <28-Aug-2013 16:35:18 tdo>
# $Id: $
#
# Copyright (C) 2013 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#


open F, shift or die "Print give me list";
my %h;

while (<F>) {
  chomp;
  my @ar=split(/\t/);
  $h{$ar[0]}{$ar[1]}=1;
}
close(F);


my %R;
my %names;

while (<STDIN>) {
  chomp;
  
  if (/^#/) {
	if (/#CHROM/) {
	  my @ar=split(/\t/);
	  for my $i (9..368) {
		$names{$i}=$ar[$i]
	  }
	}
  }
  else {
	my @ar=split(/\t/);
	if (defined($h{$ar[0]}{$ar[1]})) {
	  my $pos=$ar[0].".".$ar[1];
	  for my $i (9..368) {
		my ($dp,$old,$new)= $ar[$i] =~ /\S+:(\d+):(\d+),(\d+)/;
#		print "$names{$i}  $dp,$old,$new \t";
		
		if ($dp<10) {
		  $R{$pos}{$names{$i}}=-1
		}
		elsif ($old>=(3*$new)) {
		  $R{$pos}{$names{$i}}=0
		}
		elsif ((3*$old)<$new) {
		  $R{$pos}{$names{$i}}=1
		}
		else {
		  $R{$pos}{$names{$i}}=2
		}
#		print " $R{$pos}{$names{$i}}\n";
		
	  }
	}
  }
}

#print "Genomes";
foreach my $gx (sort {$a <=> $b}  keys  %names) {
  #print "\t$names{$gx}";
}
print "\n";

foreach my $gx (sort {$a <=> $b}  keys %names) {
  my $gX=$names{$gx};  
  #print "$gX\t";
  foreach my $gy (sort {$a <=> $b} keys %names) {
	my $gY=$names{$gy};
	#print "$gX\t$gY\t$gx\t$gy\n";
	if ($gx<$gy){
	  ### get now all, count equal, ignore -1
	  my $countOK=0;
	  my $same=0;	
	  foreach my $pos (keys %R) {
		#	  print "$pos $gX $R{$pos}{$gX}\n";
		
		if ($R{$pos}{$gX} != -1 && $R{$pos}{$gY} != -1) {
		  $countOK++;
		  if ($R{$pos}{$gX} == $R{$pos}{$gY}){
			$same++;			
		  }		  
		  
		}
	  } # end forech SNP's
	  $res[$gx][$gx]=1;
	  $res[$gy][$gy]=1;
	  #	print "$same \t $countOK\n";
	  $res[$gy][$gx]=sprintf("%.2f", (100* $same / $countOK));
	  $res[$gx][$gy]=sprintf("%.2f", (100* $same / $countOK));
	}  
  }
  #print "\n";
  
}

foreach my $gx (sort {$a <=> $b}  keys %names) {
  print "\t$names{$gx}"
}
print "\n";

foreach my $gx (sort {$a <=> $b}  keys %names) {
  print $names{$gx};
  foreach my $gy (sort {$a <=> $b}keys %names) {
	print "\t"; 
	print $res[$gx][$gy];
  }
  print "\n";
}
