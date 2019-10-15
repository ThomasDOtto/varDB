use strict;

### looks for overlap in gff...


my %h;
my %length;
my %name;
while (<STDIN>){
	chomp;
	my @ar=split(/\s+/);
	### get the name
	foreach (@ar){
	  my ($chr,$from,$to) = $_ =~
		/^(\S+):(\d+)-(\d+)$/;
#	  print "$chr $from $to\n";
	  
#$_ =~ /gene_id \"(\S+)\"/;
#	print $ar[8]."- $id\n";
#	$name{$id}=$ar[8];
		for my $i ($from..$to){
			$h{$chr}{$i}=1;
		}	
	}

}
#lose(F);

my $n;
my $last;

foreach my $chr (keys %h){
  my $start='';
  
  foreach my $pos (sort {$a <=> $b } keys  %{$h{$chr}}) {
	if ( $start eq '') {
	  $start=$pos;
	  $last =$pos
	}
	elsif ($pos > ($last+2)) {
#	  print "$pos : $start : $last \n";
	  
	  $n.= " $chr:$start-$last";
	  $start=$pos; $last=$pos;
	  
	}
	else {
	  $last=$pos;
	}
  }
  
  if ($start ne '') {
	  $n.= " $chr:$start-$last";
	}
}

print $n."\n";

