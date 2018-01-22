#!/bin/bash

samples=$1

PROJDIR=$PWD

cd $PROJDIR

while read line
do
    if [[ ! $line =~ "SampleID" ]] && [[ ! $line =~ "#" ]]
    then
	### sample details
	sampleId=`echo $line | cut -d, -f1`
	reference=`echo $line | cut -d, -f2`
	collectionPathUri=`echo $line | cut -d, -f6`

	### preview run info
	echo ""
	echo "sampleId=$sampleId"
	echo "reference=$reference"
	echo "collectionPathUri=$collectionPathUri"

	### output directory
	rundir=`printf "%s/samples/%05i" $PROJDIR $sampleId`
	rm -rf "$rundir"
	mkdir -p "$rundir"

	### run demultiplexing script
	qsub \
	    -v root="$PROJDIR",rundir="$rundir",rname="$reference",collectionPathUri="$collectionPathUri",instrument="RSII" \
	    -N "demux$sampleId" \
	    -o "$rundir"/workflow.log \
	    -j yes \
	    "$PROJDIR"/scripts/workflow-demux.sh
    fi
done < $samples
