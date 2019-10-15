#!/bin/bash

### will annotate the 
#set -e
genome=$1;
readroot=$2
name=$3

tmp=$$;
olddir=$PWD
dir="/tmp/tdo/dir.tdo.$tmp/"
#dir="local"
echo "Old dir $olddir   New Dir $dir" 

mkdir -p $dir
cp $genome $dir
cd $dir
cp -s $olddir/$readroot*q .
AUGUSTUS_CONFIG_PATH=/nfs/pathogen003/tdo/Tools/Augustus/config; export AUGUSTUS_CONFIG_PATH
### Getting some annotation. So far based on augustus
augustus --species=3D7_genefamilies $genome  > out.augustus.txt
### rename the output to the correct ID
cat out.augustus.txt | perl -e 'my $n=shift; while(<STDIN>){s/g(\d+)/$n.g$1/g; print }' $name > out2.augustus.txt

#getAnnoFasta.pl out2.augustus.txt; mv out.augustus.aa Genes.$name.aa.fasta
### to embl
~/Bin/annotation.wrapperAnnotateAugustus.sh $genome /nfs/pathogen003/tdo/Pfalciparum/3D7/Reference/Oct2011/Pf3D7_v3.aa.fasta out2.augustus.txt

mkdir Seq
cd Seq
SeperateSequences.pl ../$genome  &> /dev/null
cd ..
mkdir EMBL
for x in ` grep '>' $genome | sed 's/>//g' ` ; do perl ~/Bin/ratt/main.ratt.pl  doEMBL EMBL/$x Augustos/$x.embl Seq/$x; done &> Annotation.txt
#psu_union.pl EMBL/*embl > $name.union.embl
rm -rf Seq/

### get a union file and all the genoems with VAR annotation
list=$(ls -S EMBL/* | awk '{s=s" "$1} END{print s}')
psu_union.pl $list > Union.$name.embl
rm -rf EMBL/;
rm -rf Augustos*/
#perl /nfs/users/nfs_t/tdo/Bin/EMBL2aa_moreDescription.pl Union.$name.embl > $name.all.aa.fasta

touch All.$name.aa.fasta All.$name.na.fasta

~tdo/Bin/little.allEMBL2xx.sh $name
mv All.$name.aa.fasta $name.all.aa.fasta
mv All.$name.na.fasta $name.all.nt.fasta
#perl /nfs/users/nfs_t/tdo/Bin/EMBL2CDS.Reichenowi.previousID.pl Union.$name.embl $name.all.nt.fasta 0
grep -i PfEMP1 $name.all.aa.fasta | awk '{ print  $1 }' | sed 's/>//g' >  id.withPfEMP1.annotation.txt

### get just the genes annotated with Pfemp1
fasta_extractor.pl -i $name.all.aa.fasta -f id.withPfEMP1.annotation.txt -s $name.VAR.aa.fasta &> /dev/null
fasta_extractor.pl -i $name.all.nt.fasta -f id.withPfEMP1.annotation.txt -s $name.VAR.nt.fasta &> /dev/null


### get sis data
fasta2singleLine.pl $name.VAR.aa.fasta | perl -nle 'if (/>/){print } elsif(/(DI\wD[I|V]\wR.{10,160}DY\wPQ\w{2}R)/) { print $1} else { print "NA"};' | grep -B 1 DI | grep -v "^\-\-" > dbl.$name.fasta 
perl ~/Bin/DSID_cys.pl dbl.$name.fasta cys.$name.txt  &> /dev/null
cut -f 3 cys.$name.txt | sort | grep -v cys | uniq -c  | perl -nle 's/\s+/\t/g; print '> cys.$name.summary.txt

### differential expression
# it will be jsut done on the ..nt...
# map reads with smalt
SMALT_PARAMETER=" -n 1 -y 0.7 " ; export SMALT_PARAMETER
genome1k=$name.VAR.nt.0.5kb.fasta
genome=$name.VAR.nt.fasta
/nfs/users/nfs_t/tdo/Bin//Countcountigsize.pl $genome 500 $genome1k

~/Bin/little.smalt.bam.sh $genome1k 20 8 $readroot\_1.fastq $readroot\_2.fastq Res.Mapped 1000
/software/java/bin/java -XX:MaxPermSize=512m -Xmx2000m -XX:ParallelGCThreads=1 -XX:-UseParallelGC -XX:+UseSerialGC -jar $ICORN2_HOME/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=Res.Mapped.bam O=Res.Mapped.markdup.bam

## We will count the reads per var gene longer than 1kb
samtools view -F 1028 Res.Mapped.bam | cut -f 3 | uniq -c | perl -nle '/^\s*(\d+)\s+(\S+)/; print "$2\t$1" ' > ReadCount.0.5k.txt
amount=$(awk '{ s+=$2} END { print s}' ReadCount.0.5k.txt )
samtools faidx $genome1k
cat ReadCount.0.5k.txt | perl ~tdo/Bin/vargene.rawReadCountNormalization.pl $genome1k.fai  $amount > ReadCount.$name.txt
#~/Bin/fasta2gtf.pl $genome > $genome.gtf
#awk '$2=="AUGUSTUS"' out2.augustus.txt | egrep "transcript|exon" | sed 's/CDS/exon/g' | perl -e 'while(<STDIN>){ if (/exon/){@ar=split(/\t/); $CDS{$ar[0]}++; chomp($_); print "$_ exon_number \"$CDS{$ar[0]}\"\n"  }}' > ForCufflinks.gtf

#bam=Res.Mapped.bam
#cufflinks  --no-update-check -b $genome -q  -G $genome.gtf -o cuff.coverage $bam


### blast of the UPS sequences / interal regions.
formatdb -p F -i $genome
blastall -m 8 -F F -p blastn -e 1e-2 -d $genome -i /nfs/pathdata2/Plasmodium/falciparum/3D7/Reference/ups.sequences.txt -o comp.ups.$name.blast
blastall -m 8 -p blastn -e 1e-80 -d $genome -i /nfs/pathdata2/Plasmodium/falciparum/3D7/Reference/PfV2.internalVAR.fasta -o comp.internal.$name.blast



cp Read* cys*  *VAR*a *.all.*.fasta  Union*embl   $olddir

rm out*txt
rm Union*embl.??.fasta *fai
rm $genome* formatdb.log  id.withPfEMP1.annotation.txt
rm *fastq 
cp Res.Mapped.bam $olddir
zip -r $olddir/Annotation.zip *

echo "Old dir $olddir   New Dir $dir" 

cd $olddir
rm -rf $dir/
echo "done."
 

