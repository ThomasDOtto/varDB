

### circos renameing

my %h;
$h{ATSA}="black";
$h{ATSB}="black";
$h{CIDRa}="greys-9-seq-9";
$h{CIDRb}="greys-9-seq-8";
$h{CIDRd}="greys-9-seq-7";
$h{CIDRg}="greys-9-seq-6";
$h{DBLa}="greys-9-seq-5";
$h{DBLb}="greys-9-seq-6";
$h{DBLd}="greys-9-seq-7";
$h{DBLe}="greys-9-seq-8";
$h{DBLg}="greys-9-seq-9";
$h{DBLz}="greys-9-seq-5";
$h{NTSA}="black";
$h{NTSB}="black";

while (<>){
my @ar=split(/\t/);
my ($id) = $ar[8] =~ /Domain\.(\D{4,5})\d+/;
my $col = $h{$id};	
	print "$ar[0] $ar[3] $ar[4] fill_color=$col\n";
}