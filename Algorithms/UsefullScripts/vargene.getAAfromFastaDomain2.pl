use strict;

my $ref_h=loadFasta(shift);
my %h;

while (<STDIN>){
	chomp;
	my @ar=split(/\t/);
	$h{$ar[1]}.=">$ar[0].$ar[4]\n";
	if (!defined($$ref_h{$ar[0]})){
		print "$ar[0] - $_\n";	
	}	
	$h{$ar[1]}.=substr($$ref_h{$ar[0]},($ar[2]-1),($ar[3]-$ar[2]+1))."\n";
	
}

### write now the files…
foreach my $k (keys %h){
	open F, " > Domain.$k.fasta " or die "Problems\n";
	print F $h{$k};
	close(F);	
}

####################
### loadfasta
####################
sub loadFasta{
  my $file = shift;

  my $name;
  my %h;
  
  open (F,$file) or
        die "Couldn't open Sequence file $file: $!\n";
  
  while (<F>) {
        if (/^>(\S+)/) {
          $name=$1;
        }
        else {
          chomp;
          $h{$name}.=$_
        }
  }
  return \%h
}

