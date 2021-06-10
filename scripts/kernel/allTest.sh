#!/bin/bash

ROOT=/home/xjj_stencil
# ROOT=`pwd`
HEADINSERTER=${ROOT}/myHeadInserter
TRANSLATOR=${ROOT}/myTranslator

PHYSIS_HOME=/home/physis_home/install
PHYSIS_INCLUDE=${PHYSIS_HOME}/include

LDFLAGS="-I${PHYSIS_INCLUDE}"
RFLAGS=-DPHYSIS_USER

BENCH_HOME=/home/benchmarks-generate/kernelBench
LDFLAGS+=" -I${BENCH_HOME}"
BENCH_LIST="mintBench patusBench"



TIMES=10
RESULT=result-physis.txt

rm -rf ${RESULT}

############################################

for BENCH_NAME in ${BENCH_LIST}
do
	cd ${BENCH_NAME}

	SRC_LIST=`ls *.c`
	for SRC in ${SRC_LIST}
	do
	    EXE=${SRC%*.c}
	    echo running ${EXE}
	    for ((i = 1; i <= TIMES; i++))
	    do
	        echo round ${i}
	        ./${EXE}
	    done
	    echo
	done
	cd ../
done

