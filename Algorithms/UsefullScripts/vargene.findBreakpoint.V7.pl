#! /usr/bin/perl -w
#
# File: vargene.blastInsertSites.pl
# Time-stamp: <19-Jan-2015 16:03:43 tdo>
# $Id: $
#
# Copyright (C) 2012 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#

use strict;
use warnings;
use Carp;
use Data::Dumper;

my $type;
my $hitsFile=shift;
my $geneFile=shift;
my $BORDER=200;
my $DIST_HITS=750; 
my $showLength=15;

open F, $geneFile or die "Problem give me file: $1\n";
my %g;
my %h;
my %h2;
my %hits;

my $n;
while (<F>) {
  chomp;
  if (/>(\S+)/) {
	$n=$1
  }
  else {
	$g{$n}.=$_
  }
}
close(F);

my ($ref_hits2,$ref_hits2Sum,$ref_hitsID)=loadHits($hitsFile);

print "#loaded files\n";
while (<STDIN>) {
# PA0050_VAR_12 PH0026_VAR_3 99.88 830 1 0 6704 7533 4414 5243 0.0
# 1637 7751 5414
  chomp;
  my $rev=0;
  my ($q,$s,$id,$ov,$one,$two,$qStart,$qEnd,$sStart,$sEnd,
	  $score,$value,$qLength,$sLength)=split(/\s+/);
  if ($sEnd<$sStart) {
	my $tmp=$sEnd;$sEnd=$sStart;$sStart=$tmp;	
	$rev=1;
	warn "should not happen… $_ \n";
  }
  $hits{$s}++;
#  print Dumperxs $$ref_hits{$s};
  ### alignment is on right hand site
  if ($qStart>$BORDER && ($qStart+$BORDER)<$qLength &&
	  $sStart>$BORDER && ($sStart+$BORDER)<$sLength){
	# round it by 10 bases

	##check sequence for N's
  #      my $seq2= substr($g{$s},($sStart-1),($showLength));
  #      my $seq1= substr($g{$q},($qStart-1),($showLength));
	my $seq1= substr($g{$s},($sStart-$showLength),(2*$showLength));
	my $seq2= substr($g{$q},($qStart-$showLength),(2*$showLength));
	
	##check sequence for N's
	if ( !($seq1=~/N/i || $seq2=~/N/i) ) {

	    my $leftQuery=0;
	    my $leftSubject=0;
	    if (defined($$ref_hits2Sum{$s}{($sStart-1)}{$q})){
		$leftSubject=$$ref_hits2Sum{$s}{($sStart-1)}{$q};
	    }
	    if (defined($$ref_hits2Sum{$q}{($qStart-1)}{$s})){
		$leftQuery=$$ref_hits2Sum{$q}{($qStart-1)}{$s};
	    }

	    my $leftSumQuery=0;
	    my $leftSumSubject=0;
	    my ($genomeQuery) = $q =~ /^(\S+-C\S*)\.g/;
	    if (!defined($genomeQuery)){
		($genomeQuery) = $q =~ /^(\S+)_\d+/;
	    }
	    my ($genomeSubject) = $s =~ /^(\S+-C\S*)\.g/;
	    if (!defined($genomeSubject)){
		($genomeSubject) = $s =~ /^(\S+)_\d+/;
            }
	    if (!defined($genomeQuery)){
		die "$_\n";
	    }
	    if (defined($$ref_hits2Sum{$q}{($qStart-1)})){
		$leftSumQuery   = scalar (grep { /$genomeQuery/ } keys %{$$ref_hits2Sum{$q}{($qStart-1)}});
	    }
	    if (defined($$ref_hits2Sum{$s}{($sStart-1)})){
		$leftSumSubject = scalar (grep { /$genomeSubject/ } keys %{$$ref_hits2Sum{$s}{($sStart-1)}});
	    }

	    if ($$ref_hits2Sum{$s}{$sStart}{$q}>=1  && 
		( (defined($$ref_hits2Sum{$s}{$sStart}{$q} ) && $$ref_hits2Sum{$s}{$sStart}{$q} > 1) 
		  || (defined($$ref_hits2Sum{$q}{$qStart}{$s}) && $$ref_hits2Sum{$q}{$qStart}{$s} >=2)
		)
		){
		
		$type="Repeat\t";
	    }
	    elsif( 
		 ( ($leftSumQuery-$leftQuery)>=1)
		||
		   (($leftSumSubject-$leftSubject)>=1   )

#scalar(keys $$ref_hits2Sum{$s}{$sStart})>1 &&
#                   ( scalar(keys $$ref_hits2Sum{$s}{($sStart-1)})>=1
#                     || scalar(keys $$ref_hits2Sum{$q}{($qStart-1)})>=1)
                   )
            {
                $type="DuplicatedGene\t";
            }
	    # contained
	    elsif ($$ref_hits2Sum{$s}{$sStart}{$q}>1 &&  defined($leftQuery) &&  ($sStart-$leftQuery)==1  ){
		$type="Contained\t";
	    }
	    # indel subject
	    elsif($$ref_hits2Sum{$s}{$sStart}{$q}==1 && $leftQuery==0 &&  $leftSubject>=1 ){
		$type="Insert\t";
	    }
	    elsif($$ref_hits2Sum{$s}{$sStart}{$q}==1 
		  && $leftQuery>=1 &&  $leftSubject==0 ){
		$type="Deletion\t";
	    }
	    elsif($$ref_hits2Sum{$s}{$sStart}{$q}==1
		  && $leftQuery==0
		  && $leftSubject==0 
## chedk one hits
		  && defined($$ref_hitsID{$s.":::".($sStart-$DIST_HITS)}{$q})
            ){
            $type="DOUBLE ".$$ref_hitsID{$s.":::".($sStart-$DIST_HITS)}{$q}."\t";
        }
       elsif($$ref_hits2Sum{$s}{$sStart}{$q}==1 
	      && $leftQuery==0 
	      && $leftSubject==0
	    ){
	    $type="OK\t";
        }
	    else {
		$type="IGNORE\t";
	    }
	    print $type;	  
	    
	    my $res= "$id\t$ov\t$qStart-$qEnd\t$sStart-$sEnd\t$qLength\t$sLength - LEFT";
	    print "$seq1\t$q\t$qStart\t$s\t$sStart\t$res\n";
#	  print "$seq1 -  $res - A\n";
#	  print "$seq2 - $res - A\n";
	    my $pos=$sStart;#(12* (int ($sStart/12)));
	    $h{$s}{$pos}++;
	    $h2{$s}{$pos}.= "$q\t$id\t$ov\t$qStart\t$qEnd\t$sStart\t$sEnd\n";
	}
	
  }
  
  ### alignment is on left hand site, (end of alignem
  if ($qEnd>$BORDER && ($qEnd+$BORDER)<$qLength &&
	  $sEnd>$BORDER && ($sEnd+$BORDER)<$sLength){
      # round it by 10 bases
      my $seq2= substr($g{$s},($sEnd-$showLength-1),(2*$showLength));
	my $seq1= substr($g{$q},($qEnd-$showLength-1),(2*$showLength));
	
  	if ( !($seq1=~/N/i || $seq2=~/N/i) ) {


	    my $rightQuery=0;
	    my $rightSubject=0;
	    if (defined($$ref_hits2Sum{$s}{($sEnd+1)}{$q})){
		$rightSubject=$$ref_hits2Sum{$s}{($sEnd+1)}{$q};
	    }
	    if (defined($$ref_hits2Sum{$q}{($qEnd+1)}{$s})){
		$rightQuery=$$ref_hits2Sum{$q}{($qEnd+1)}{$s};
	    }


	    my $rightSumQuery=0;
            my $rightSumSubject=0;
            my ($genomeQuery) = $q =~ /^(\S+-C\S*)\.g/;
            if (!defined($genomeQuery)){
                ($genomeQuery) = $q =~ /^(\S+)_\d+/;
            }
            my ($genomeSubject) = $s =~ /^(\S+-C\S*)\.g/;
            if (!defined($genomeSubject)){
                ($genomeSubject) = $s =~ /^(\S+)_\d+/;
            }

            if (defined($$ref_hits2Sum{$q}{($qEnd+1)})){
                $rightSumQuery   = scalar (grep { /$genomeQuery/ } keys %{$$ref_hits2Sum{$q}{($qEnd-1)}});
            }
            if (defined($$ref_hits2Sum{$s}{($sEnd-1)})){
                $rightSumSubject = scalar (grep { /$genomeSubject/ } keys %{$$ref_hits2Sum{$s}{($sEnd-1)}});
            }

	    if (!defined($genomeQuery)){
		die "genomeQuery $_\n";
	    }
	    ## Repaeat
  	   if ($$ref_hits2Sum{$s}{$sEnd}{$q} >= 1 
	       && 
	       ( (defined($$ref_hits2Sum{$s}{$sEnd}{$q}) && $$ref_hits2Sum{$s}{$sEnd}{$q} > 1) 
		 || ( defined($$ref_hits2Sum{$q}{$qEnd}{$s}) && $$ref_hits2Sum{$q}{$qEnd}{$s} >=2) )
	       ) {
#defined($$ref_hits2Sum{$s}{($sEnd+1)}{$q.":::".($qEnd+3)}) && $$ref_hits2Sum{$s}{($sEnd+1)}{$q.":::".($qEnd+3)}>1) { 
	       $type="Repeat\t";
	   }
	    elsif( 
		 ( ($rightSumQuery-$rightQuery)>=1)
		||
		   (($rightSumSubject-$rightSubject)>=1   )
		)

	    {
		$type="DuplicatedGene\t";
	    }
	   # contained	 
	   elsif ($$ref_hits2Sum{$s}{$sEnd}{$q}>1 &&  defined($$ref_hits2Sum{$s}{($sEnd+1)}{$q}) 
		  &&   $$ref_hits2Sum{$s}{($sEnd+1)}{$q}==1){
	       $type="Contained\t";
	   }
	   # indel subject	
	    # indel subject
	    elsif($$ref_hits2Sum{$s}{$sEnd}{$q}==1 && $rightQuery==0 &&  $rightSubject>=1 ){
		$type="Insert\t";
	    }
	    elsif($$ref_hits2Sum{$s}{$sEnd}{$q}==1 
		  && $rightQuery>=1 &&  $rightSubject==0 ){
		$type="Deletion\t";
	}

#	   elsif($$ref_hits2Sum{$s}{$sEnd}{$q}>=1 
#		 && defined($$ref_hits2Sum{$s}{($sEnd+1)}{$q})  
#		 && !defined($$ref_hits2Sum{$q}{($qEnd+1)}{$s} )
#	       ){
#	       $type="Insert\t";
#	   }
#	   elsif($$ref_hits2Sum{$s}{$sEnd}{$q}>=1 
#		 && defined($$ref_hits2{$q.":::".($sEnd+1)}{$s} 
#			    && ($qEnd - $$ref_hits2{$q.":::".($sEnd+1)}{$s})<250) ){
#	       $type="Deletion\t";
#	   }
	    elsif($$ref_hits2Sum{$s}{$sEnd}{$q}>=1
		  && $rightSubject==0
		  && $rightQuery==0
		  && defined($$ref_hitsID{$s.":::".($sEnd+$DIST_HITS)}{$q})
		){
		$type="DOUBLE ".$$ref_hitsID{$s.":::".($sEnd+$DIST_HITS)}{$q}."\t";
	    }
	    elsif($$ref_hits2Sum{$s}{$sEnd}{$q}==1 
		 && $rightSubject==0
		 && $rightQuery==0
	       ){
	       $type="OK\t";
	   }

	   else {
	       $type="IGNORE\t";
	   }
	   print $type;
	   my $res= "$id\t$ov\t$qStart-$qEnd\t$sStart-$sEnd\t$qLength\t$sLength - RIGHT";
	   print "$seq1\t$q\t$qEnd\t$s\t$sEnd\t$res\n";
#	  print "$seq1 - $res - E\n";
	   
#	  print "$seq2 - $res - E\n";
	   my $pos=$sEnd;#(12* (int ($sEnd/12)));
	   $h{$s}{$pos}++;
	   $h2{$s}{$pos}.= "$q\t$id\t$ov\t$qStart\t$qEnd\t$sStart\t$sEnd\n";
	}
  }
}

if (0){
foreach my $s (keys %h) {
#xs  print ">$s: ($hits2Sum{$s} hits)\n";
  
  foreach my $pos (sort {$a<=> $b} keys %{$h{$s}}) {
	### start is always in one straind
	my $posN=($pos - (($pos-1)%3)-1);
	
	my $ref = sixFrameTranslate(substr($g{$s},($posN-30),60));
	
	print "$pos - $posN -  $h{$s}{$pos}\t".substr($g{$s},($pos-10),20)."\t".$$ref[0];
#	foreach (@$ref) {
#	  print "\t$_";
#  
#	}
	print "\n";
	
  }
}
}
sub revCompDNA
{
        croak "Usage: revCompDNA string\n" unless @_ == 1;
        my $seq = shift;
        $seq= uc( $seq );
        $seq=reverse( $seq );
        
        $seq =~ tr/ACGTacgt/TGCAtgca/;
        return $seq;
}

sub loadHits
{
my %hitsID;
my %hits2;
my %hits2Sum;

open F, shift or die "Hit file not found... $!\n";

while (<F>){
   my @ar=split(/\t/);
   if ($ar[0] ne $ar[1] && $ar[2]>=90 &&
	$ar[3]>=250){

#PF3D7_0100100	PM0087-C.g463	99.82	1110	2	0	3706	4815	3460	4569	0.0	2184	6492	4743
       for my $i ($ar[8] .. $ar[9]){
	   my $c=($ar[6]+$i-$ar[8]);
	   
	   $hitsID{$ar[1].":::".$i}{$ar[0]}=$ar[2];
	   $hits2{$ar[1].":::".$i}{$ar[0]}=$c;
	   $hits2Sum{$ar[1]}{$i}{$ar[0]}++;
       } 
       ## if not 100%, does not add up
   

   }
   
}
return (\%hits2,\%hits2Sum,\%hitsID);;
}

sub sixFrameTranslate
{
        croak "Usage: sixFrameTranslate string\n" unless @_ == 1;
        my $forward = $_[ 0 ];
        my $reverse = revCompDNA( $_[ 0 ] );
        my @translations;
        
        #forward 3 frame translation
        for( my $i = 0; $i < 1; $i ++ )
        {
                my $aa = '';
                my $j = $i;
                while()
                {
                        if( length( $forward ) < $j + 3 )
                        #if( ! defined substr( $forward, $i, 2 ) )
                        {
                                last;
                        }
                        else
                        {
                                $aa .= codon2aa( substr( $forward, $j, 3 ) );
                        }
                        $j += 3;
                }
                push( @translations, $aa );
        }
        
        #reverse 3 frame translation
        for( my $i = 0; $i < 0; $i ++ )
        {
                my $aa = '';
                my $j = $i;
                while()
                {
                        if( length( $reverse ) < $j + 3 )
                        #if( ! defined $reverse[ $i + 2 ] )
                        {
                                last;
                        }
                        else
                        {
                                $aa .= codon2aa( substr( $reverse, $j, 3 ) );
                        }
                        $j += 3;
                }
                push( @translations, $aa );
        }
        
        return \@translations;
}

sub codon2aa
{
        croak "Usage: codon2aa codon\n" unless @_ == 1;
        my( $codon ) = @_;
                
		my %genetic_code = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
);
                
        $codon = uc $codon;
                
        if( exists $genetic_code{ $codon } )
        {
                return $genetic_code{ $codon };
        }
        else
        {
                return '*';
        }
}

									
