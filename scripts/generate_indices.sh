#!/bin/bash

organism=$1
taxid=$2
gtf=$3
fasta=$4
appendix=$5

outdir='/home/data/indices'

DEXSeqScript='/home/scripts/DEXSeq/dexseq_prepare_annotation.py'
RIndices='/home/scripts/generate_R_index.R'
SalmonScript='/home/software/gffread/gffread-0.11.5.Linux_x86_64/gffread'
KallistoScript='/home/software/kallisto_linux-v0.45.0/kallisto'
STARScript='/home/software/STAR/bin/Linux_x86_64_static/STAR'
HISATScript='/home/scripts/hisat2_index.sh'

echo 'Organism'$'\t'$organism
echo 'taxid'$'\t'$taxid
echo 'outDir'$'\t'$outdir
echo 'appendix'$'\t'$appendix
echo 'gtf'$'\t'$gtf
echo 'fasta'$'\t'$fasta

#DEXSeq
echo $'\n'"generating DEXSeq index..."
mkdir -p $outdir/DEXSeq/$appendix/
/usr/bin/python3 $DEXSeqScript --aggregate=no $gtf $outdir/DEXSeq/$appendix/annot.noaggregate.gtf
echo "second DEXSeq gtf..."
/usr/bin/python3 $DEXSeqScript $gtf $outdir/DEXSeq/$appendix/annot.gtf

#R
echo $'\n'"generating R index..."
mkdir -p $outdir/R/$appendix/
$RIndices --gtf $gtf --outdir $outdir/R/$appendix/ --organism $organism --taxonomyId $taxid


#SALMON
echo $'\n'"generating salmon index..."
mkdir -p $outdir/salmon/$appendix/
echo "gffread -w $outdir/salmon/$appendix/cdna.fa -g $fasta $gtf"
$SalmonScript -w $outdir/salmon/$appendix/cdna.fa -g $fasta $gtf


#KALLISTO
echo $'\n'"generating kallisto index..."
mkdir -p $outdir/kallisto/$appendix/
echo "kallisto index --index $outdir/kallisto/$appendix/INDEX $outdir/salmon/$appendix/cdna.fa"
$KallistoScript index --index $outdir/kallisto/$appendix/INDEX $outdir/salmon/$appendix/cdna.fa

#HISAT2
echo $'\n'"generating HISAT2 index..."
mkdir -p $outdir/hisat2/$appendix/
echo "hisat2_index.sh $gtf $fasta $outdir/hisat2/$appendix/ $outdir/hisat2/$appendix/INDEX"
$HISATScript $gtf $fasta $outdir/hisat2/$appendix/ $outdir/hisat2/$appendix/INDEX

#STAR
echo $'\n'"generating STAR index..."
overhang=99

mkdir -p $outdir/STAR/$appendix/INDEX/$overhang
echo "$STARScript --runMode genomeGenerate --runThreadN 8 \
    --genomeDir $outdir/STAR/$appendix/INDEX/ --genomeFastaFiles $fasta \
    --sjdbGTFfile $gtf --sjdbOverhang $overhang"

$STARScript --runMode genomeGenerate --runThreadN 8 \
    --genomeDir $outdir/STAR/$appendix/INDEX/ --genomeFastaFiles $fasta \
    --sjdbGTFfile $gtf --sjdbOverhang $overhang
