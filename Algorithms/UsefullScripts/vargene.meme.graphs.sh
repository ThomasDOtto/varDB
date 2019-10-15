

res=$1
num=$2
grep -A 10000000 "SEQUENCE.*COMBINED.*P-VALUE" meme.txt | perl -e 'while (<>){ chomp; $n.=$_; if (/\\$/){ } else { $n=~ s/\\\s+//g; print $n."\n"; $n=""}};' | egrep "^DBL|^ATS|^NTS|^CID|^P" | perl ~tdo/Bin/avian.parseMeme.pl $num | perl ~tdo/Bin/vargene.renameMatrix.DBLstuff.v2.pl ../List.Rename.Domain.txt > Matrix.$res.txt
cut -f 1 Matrix.$res.txt  | egrep -v Name | perl -nle 'if (/(CIDR\S\d\.\d+)/) { print $1 } elsif (/(CIDR.{1,6})\./) { print $1 } elsif (/(DBL\S\d\.\d+)/) { print $1} elsif (/(DBL.{0,6})\./) { print $1} elsif (/(.TS.{0,3}\d+)/) { print $1} else { print "UNDEF"} ' | perl ~tdo/Bin/heatmap.List2Matrix.feature_selfIndex.pl > DBL.$res.txt
R CMD BATCH "--args Matrix.$res.txt DBL.$res.txt DBLSub.$res " ~tdo/Bin/DoHeatmap.Species.right.R
