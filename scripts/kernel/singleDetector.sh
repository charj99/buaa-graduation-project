#!/bin/bash

ROOT=/home/xjj_stencil
# ROOT=`pwd`
HEADINSERTER=${ROOT}/myHeadInserter
TRANSLATOR=${ROOT}/myTranslator

PHYSIS_HOME=/home/physis_home/install
PHYSIS_INCLUDE=${PHYSIS_HOME}/include

LDFLAGS="-I${PHYSIS_INCLUDE}"
RFLAGS=-DPHYSIS_USER

if [ $# -lt 2 ]
then
    echo "please provide input file and benchmark name"
    exit
fi

BENCH_NAME=$2
BENCH_HOME=/home/benchmarks-generate/kernelBench/${BENCH_NAME}
LDFLAGS+=" -I${BENCH_HOME}/.."

if [ ! -d ${BENCH_NAME} ]
then
    mkdir ${BENCH_NAME}
fi

INFO_NAME=${1%*.c}.info

# use headInserter and myTranslator to translate C file
${HEADINSERTER} ${BENCH_HOME}/$1 ${LDFLAGS} ${RFLAGS}
mv ./rose_$1 ./$1
${TRANSLATOR} ./$1 ${LDFLAGS} ${RFLAGS} 2>&1 | tee ./${INFO_NAME}

# save useful files and remove useless files
mv ./rose_$1 ./${BENCH_NAME}
mv ./${INFO_NAME} ./${BENCH_NAME}
rm -rf $1
