#!/bin/bash
ROOT=/home/xjj_stencil
# ROOT=`pwd`
HEADINSERTER=${ROOT}/myHeadInserter
TRANSLATOR=${ROOT}/myTranslator

PHYSIS_HOME=/home/physis_home/install
PHYSIS_INCLUDE=${PHYSIS_HOME}/include
PHYSIS_LIB=${PHYSIS_HOME}/lib
PHYSISC=${PHYSIS_HOME}/bin/physisc

LDFLAGS="-I${PHYSIS_INCLUDE}"
RFLAGS=-DPHYSIS_USER

BENCH_HOME=/home/benchmarks-generate/kernelBench
LDFLAGS+=" -I${BENCH_HOME}"
BENCH_LIST="mintBench patusBench"

CUDA_ROOT=/usr/local/cuda
CUDA_BIN=${CUDA_ROOT}/bin
NVCC=${CUDA_BIN}/nvcc

PFLAGS="-I${PHYSIS_INCLUDE} -I${BENCH_HOME}"


for BENCH_NAME in ${BENCH_LIST}
do
	if [ ! -d ${BENCH_NAME} ]
	then
	    mkdir ${BENCH_NAME}
	fi
	for DIR_FILE in `ls ${BENCH_HOME}/${BENCH_NAME}/*.c`
    do
	    FILE=${DIR_FILE##*/}
		ROSE_FILE=rose_${FILE}
		TARGET_FILE=${ROSE_FILE%*.c}
		CUDA_FILE=${TARGET_FILE}.cuda.cu
		INFO_NAME=${FILE%*.c}.info

		# use headInserter and myTranslator to translate C file
		${HEADINSERTER} ${BENCH_HOME}/${BENCH_NAME}/${FILE} ${LDFLAGS} ${RFLAGS}
		mv ./${ROSE_FILE} ./${FILE}
		${TRANSLATOR} -m ./${FILE} ${LDFLAGS} ${RFLAGS} 2>&1 | tee ./${INFO_NAME}

		# save useful files and remove useless files
		mv ./${ROSE_FILE} ./${BENCH_NAME}
		mv ./${INFO_NAME} ./${BENCH_NAME}
		rm -rf ${FILE}

		# run physis to do the translation
		${PHYSISC} ${PFLAGS} --cuda ${BENCH_NAME}/${ROSE_FILE}

		# run nvcc to compile the *.cuda.cu file
		${NVCC} ${PFLAGS} -L${PHYSIS_LIB} -lphysis_rt_cuda -lphysis_rt_mpi_cuda -lphysis_rt_mpi -lphysis_rt_ref ${CUDA_FILE} -o ${TARGET_FILE}

		mv ${CUDA_FILE} ./${BENCH_NAME}
		mv ${TARGET_FILE} ./${BENCH_NAME}
    done
done
