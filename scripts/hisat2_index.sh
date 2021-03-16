#!/bin/bash
GTF=$1
FA=$2
TMP=$3
OUTBASE=$4

HISAT="/home/software/hisat2/"

echo "extracting splice sites..."
python2 $HISAT/extract_splice_sites.py $GTF >> $TMP/tmp.ss

echo "extracting exons..."
python2 $HISAT/extract_exons.py $GTF >> $TMP/tmp.exon

echo "building index..."
python2 $HISAT/hisat2-build -p 8 --ss $TMP/tmp.ss --exon $TMP/tmp.exon $FA $OUTBASE

rm $TMP/tmp.ss
rm $TMP/tmp.exon
