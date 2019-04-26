#!/bin/bash

if [ $# -lt 1 ]
then
    echo "NEED ONE ARGUMENT: RS"
    exit
fi

RS=$1

if [ $RS != 57 -a $RS != 62 -a $RS != 67 -a $RS != 70 ]
then
    echo "ROADSETS SHOULD BE 57, 62, 67, OR 70"
    exit
fi

# NBINS=15
# for FITRANGE in 60000 50000 40000 30000 20000 10000
# do
#     CUT=4.2
#     CA=NO
#     echo $RS $NBINS $FITRANGE $CUT $CA
#     csr  $RS $NBINS $FITRANGE $CUT $CA
# done

for NBINS in 15 20 25 30
do
    for FITRANGE in 60000 50000 40000 30000 20000 10000
    do
	for CUT in 4.2 4.7 5.0 5.3
	do
	    for CA in CA NO
	    do
		echo $RS $NBINS $FITRANGE $CUT $CA
		csr $RS $NBINS $FITRANGE $CUT $CA
	    done
	done
    done
done
