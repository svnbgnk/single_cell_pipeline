#!/bin/sh -e

flexbar --reads reads.fastq.gz --target result_gz --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 0.1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.fastq result_gz.fastq`

if ! $a ; then
echo "Error testing right mode gzip fastq"
echo $a
exit 1
else
echo "Test gzip OK"
fi


flexbar --reads reads.fastq.bz2 --target result_bz2 --adapter-min-overlap 4 --adapters adapters.fasta --min-read-length 10 --adapter-error-rate 0.1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_right.fastq result_bz2.fastq`

if ! $a ; then
echo "Error testing right mode bzip2 fastq"
echo $a
exit 1
else
echo "Test bzip2 OK"
fi

echo ""

