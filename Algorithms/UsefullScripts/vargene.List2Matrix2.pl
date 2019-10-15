use strict;

my %h;
my @order;
my %lines;
my %h2;
my %c;

while (<>) {
  chomp;
  # ORTHOMCL70/f1.genes.aln:>PH0686-C.g337
  /^(\S+)\t(\S+)$/;
	$h{$2}=1;
  $h2{$1}{$2}++;
if (!defined($lines{$1})){ push @order, $1}
  $lines{$1}=1;

}

print "ID";
foreach my $k (sort { $a <=> $b} keys %h) {
  print "\t$k";
}
print "\n";


foreach my $line (@order) {
  print "$line";
  foreach my $k (sort { $a <=> $b} keys %h ) {
	if (defined($h2{$line}{$k}) ) {
	  printf ("\t%d", ($h2{$line}{$k}));
	}
	else {
	  print "\t0"
	}
  }
  print "\n";
  
}

