#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 8

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

### update shell $PATH to include required tools
module load smrtlink-5.0.0
module load samtools-1.3.1

### stop if we are dealing with unknown reference
if [ "$rname" != "b3" ] && [ "$rname" != "b4" ]; then
    echo "[ERROR] Unknown reference '$rname'"
    exit
fi

reference="$root/references/$rname/sequence/$rname.fasta"

### switch to working directory
mkdir -p "$rundir"
cd "$rundir"

### looking up sequencing data
if [ "$instrument" == "RSII" ]; then
    echo ""
    find "$collectionPathUri"/Analysis_Results -type f -name "*.bax.h5" | sort -u > input.fofn
    echo "Task 1 completed at $(date)"
else
    echo ""
    find "$collectionPathUri" -type f -name "*.subreads.bam" | sort -u > input.fofn
    echo "Task 1 completed at $(date)"
    
    echo ""
    if [ `cat input.fofn | wc -l` -eq 1 ] ; then
	ln -s `cat input.fofn` movie.subreads.bam
    else
	print "[ERROR] Cannot uniquely locate subreads.bam"
	exit
    fi
fi

### convert bax to bam
### check if FOFN contains filenames with spaces
echo ""
num=`cat input.fofn | grep ' ' | wc -l`

if [ $num -gt 0 ]
then
    ### there are spaces, use workaround
    IFS=$'\r\n' GLOBIGNORE='*' command eval 'hdf=($(<input.fofn))'
    bax2bam "${hdf[@]}" -o movie
else
    ### no spaces, proceed as usual
    bax2bam --fofn=input.fofn -o movie
fi
echo "Task 2 completed at $(date)"

### cluster reads
echo ""
if [ "$rname" == "b3" ]; then
    "$root"/bin/cluster.pl --insert_length 99 movie.subreads.bam | bzip2 - > clusters.csv.bz2
elif [ "$rname" == "b4" ]; then
    "$root"/bin/cluster.pl --insert_length 98 movie.subreads.bam | bzip2 - > clusters.csv.bz2
fi
echo "Task 3 completed at $(date)"

### split forward and reverse reads
echo ""
"$root"/bin/split.pl movie.subreads.bam clusters.csv.bz2
echo "Task 4 completed at $(date)"

for strand in fwd rev
do
    ### convert to bam
    echo ""
    samtools view -Sb subreads.${strand}.sam > subreads.${strand}.bam
    rm -f subreads.${strand}.sam
    echo "Task 5 ($strand) completed at $(date)"

    ### build ccs
    echo ""
    ccs --reportFile=subreads_ccs.${strand}.csv --logFile=subreads_ccs.${strand}.log --numThreads=8 --minPasses=1 subreads.${strand}.bam subreads_ccs.${strand}.bam
    echo "Task 6 ($strand) completed at $(date)"

    ### align reads
    echo ""
    blasr subreads_ccs.${strand}.bam "$reference" --bestn 1 --clipping soft --bam --out aligned_reads.${strand}.bam
    echo "Task 7 ($strand) completed at $(date)"

    ### extract barcodes and overhangs
    if [ "$rname" == "b3" ]; then
	"$root"/bin/extract.pl \
	    --region 1-3,1-6,7-9 \
	    --region 46-48,48-52/1,52-54 \
	    --region 91-93,94-99,97-99 \
	    aligned_reads.${strand}.bam "$reference" >fragments.${strand}.csv 2>fragments.${strand}.log
    elif [ "$rname" == "b4" ]; then
	"$root"/bin/extract.pl \
	    --region 1-3,3-10/1,10-12 \
	    --region 45-47,47-52/1,52-54 \
	    --region 87-89,89-96/1,96-98 \
	    aligned_reads.${strand}.bam "$reference" >fragments.${strand}.csv 2>fragments.${strand}.log
    fi
    echo "Task 8 ($strand) completed at $(date)"

    ### extract ZMW stats
    echo ""
    "$root"/bin/bam2csv.pl subreads_ccs.${strand}.bam zmws.${strand}.csv
    echo "Task 9 ($strand) completed at $(date)"
done

### cleanup
rm -rf movie.*
rm -rf subreads.*
echo "Task 10 completed at $(date)"

### extract results
echo ""
mkdir -p "$rundir"/results
cd "$rundir"/results
"$root"/bin/reporter.pl \
    --etype "$rname" \
    --match 5 \
    --np 5 \
    --bcout barcodes.csv \
    --ohout overhangs.csv \
    ..
echo "Task 11 completed at $(date)"

echo ""
echo "Finished on $(date)"
