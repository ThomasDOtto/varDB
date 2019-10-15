use strict;

my $ref_h=getF(shift);

while (<>){
chomp;
my @ar=split(/\t/);
if (defined($$ref_h{$ar[0]})){
$ar[0]=$$ref_h{$ar[0]}
}
print join("\t",@ar)."\n";
}

sub getF{

open F, shift or die " Problems\n";
my %h;
while (<F>){
chomp;
my @ar=split(/\t/);
if (defined($ar[1])){ $h{$ar[0]}=$ar[1];}
else {
my $n;
if (/DBL/){
($n) = $ar[0] =~ /DBL\S+([Se|P]\S+)/;
} else 
{
($n) = $ar[0] =~ /CID\S+([Se|P]\S+)/;
}
$h{$ar[0]}=$ar[1].".".$n;
}
}
return \%h;
}
