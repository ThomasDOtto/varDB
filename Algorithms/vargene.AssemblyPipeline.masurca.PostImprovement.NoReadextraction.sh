#!/bin/bash
name=$1


### Version 4: read correction
#set -e

if [ ! -d "Assembly.$name" ]; then
	echo "Direcoty does not exist";
	exit 1;
fi

cd Assembly.$name
if [ ! -f "SC.$name.fasta" ]; then
	echo "SC.$name.fasta does not exist";
	exit 1;
fi

threads=4;
SMALT_PARAMETER=" -n $threads "; export SMALT_PARAMETER
ICORN2_THREADS=$threads; export ICORN2_THREADS;


finalDIR=$(pwd)
deleteIT=0

tmp=$$;
olddir=$PWD
dir="/tmp/tdo/dir.tdo.$tmp/"
#dir="local"
echo "Old dir $olddir   New DirTEST $dir" 

mkdir -p $dir
cp SC.$name.fasta $dir
cd $dir

bam=$name.bam

### link the bam file over
echo "cp -s $olddir/$bam* ../"


# General test
insertSize=350

ln -s $olddir/rawreads_1.fastq rawreads_1.fastq
ln -s $olddir/rawreads_2.fastq rawreads_2.fastq

### 
fastq=$bam

### scaffaolding
### gapwalking Take the 
source /nfs/pathogen003/tdo/Tools/PAGIT/64bit/PAGIT/sourceme.pagit


cp SC.$name.fasta  01.assembly.fa
### Do the merging stuff of contigs, and delete potential contained...
formatdb -p F -i 01.assembly.fa
megablast -W 40 -F F -a $threads -m 8 -e 1e-80 -d 01.assembly.fa -i 01.assembly.fa | awk '$3>99 && $4>200 && $1 != $2' > comp.self1.blast
perl /nfs/users/nfs_t/tdo/Bin/helper.putlengthfasta2Blastm8.pl 01.assembly.fa 01.assembly.fa comp.self1.blast &> /dev/null



# 2a: delete contained contigs
# we want the query to be always the smaller one
awk '$3>99 && $4>500 && $13 < $14' comp.self1.blast.length | perl ~/Bin/lav.decontamination_overlapChecker.pl > 02.ListContained.txt
cat 02.ListContained.txt | awk '$4>90' | cut -f 1 > List.Contained.fofn

perl /nfs/users/nfs_t/tdo/Bin/Assembly.deleteContigs.pl List.Contained.fofn 01.assembly.fa 02.assembly.fa
### now delete the contigs from the blast
# also filter for just hits > 2kb > 99%
cat comp.self1.blast.length | awk '$3> 99.5 && $4 > 500' | perl ~tdo/Bin/pacbio.deleteEntryinBlast.pl List.Contained.fofn > Blast.merge.blast.length

# 2b: find overlaps 
### Get a bam files to do the merge
~tdo/Bin/little.smalt.bam.sh 02.assembly.fa 20 3 rawreads\_1.fastq rawreads\_2.fastq first 800 

perl ~tdo/Bin/findoverlaps_ver3.pl Blast.merge.blast.length first.bam 02.assembly.fa OUT
mv mergedseq_OUT.fasta 03.assembly.fa
rm first.bam*
###IMAGE STEPS

ln -s  03.assembly.fa Res.image.1.fasta 
### iterative image
kmer=71
mkdir IMAGE;
cd IMAGE;


cp ../Res.image.1.fasta .
$PAGIT_HOME/IMAGE/image.pl -scaffolds Res.image.1.fasta  -prefix ../rawreads -iteration 1 -all_iteration 2 -dir_prefix ite -kmer $kmer -smalt_minScore 65 &> Image.output.txt
perl /nfs/users/nfs_t/tdo/Bin/gapclosing.joinContigs.pl ite2/new.fa ite2/new.read.placed Res.image.3
rm  Image.output.txt

echo "lib2 ../rawreads_1.fastq ../rawreads_2.fastq $insertSize 0.5 FR" > lib2
perl ~tdo/bin/SSPACE-BASIC-2.0_linux-x86_64/SSPACE_Basic_v2.0.pl -l lib2  -n 31 -s Res.image.3.fasta  -x 0 -k 10 -b out.sspace > Sspace.output.txt
sed 's/|/_/g' out.sspace.final.scaffolds.fasta > ../Res.Step3.fasta
bam.correctLineLength.sh  ../Res.Step3.fasta
cd ..
if [ ! -f "Res.Step3.fasta" ] ; then
	echo "Error in Image"
	exit 127;
fi
if [ ! -s "Res.Step3.fasta" ] ; then
	echo "Image Error. File Res.image.$j.fasta is empty."
	exit 128;
fi

### clean
rm -rf IMAGE/
~tdo/Bin/assstats.sh Res.Step3.fasta > stats.IMAGED.txt

cp Res.Step3.fasta 04.assembly.fa

### icorn 2 iterations
/nfs/users/nfs_t/tdo/Bin/icorn2.serial.sh rawreads 700 Res.Step3.fasta 1 2 &> icorn2.out.txt
perl ~tdo/Bin/icorn2.collectResults.pl . > Stats.icorn2.txt
cp ICORN2.Res.Step3.fasta.3  icorned.fasta
rm -rf ICORN2*
if [ ! -f "icorned.fasta" ] ; then
	echo "Error in icorn..."
	exit 127;
fi

cp icorned.fasta 05.assembly.fa

### reapr
REAPR_AGRESSIVE_BREAKING=" -a "; export REAPR_AGRESSIVE_BREAKING
~tdo/Bin/assemblypipeline.doreapr.sh rawreads icorned.fasta SC.fixed.$name 0 $insertSize
rm -rf SC.fixed.$name/
cp  SC.fixed.$name.fasta  06.assembly.fa
### tmp stuff - moving the files over...
cp SC.fixed.$name.fasta *assembly.fa $olddir

### will rename the old reads, which is ok...

zip $olddir/resultFile.txt.zip *.txt *assembly.fa
cd $olddir
rm -rf  $dir
## go one up 
if [ -s "SC.fixed.$name.fasta" ] ; then 
	echo "Doing the annotation now... ~/Bin/vargene.AnnotationVARgenes.assembly.V0.1.sh SC.fixed.$name.fasta empty.bam $name "
	bsub.py 2 Annotate.$name bash ~tdo/Bin/vargene.AnnotationVARgenes.assembly.V0.2.cov.sh  SC.fixed.$name.fasta rawreads $name
fi
### do de novo assembly using predifined k-mer of 71

# rm *fastq *.sam


