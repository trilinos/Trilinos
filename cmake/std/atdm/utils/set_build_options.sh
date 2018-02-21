#!/bin/bash
echo Setting default compiler and build options
# Set the defaults
export COMPILER=GNU
export BUILD_TYPE=DEBUG
export USE_OPENMP=OFF
export USE_CUDA=OFF
export USE_PTHREADS=OFF
export PROJECT=PANZER

# Set the project name
if [[ $JOB_NAME == *"Trilinos"* ]]; then export PROJECT=PANZER; fi
if [[ $JOB_NAME == *"PIC"* ]];      then export PROJECT=EMPIRE-PIC; fi

# Set the compiler
if [[ $JOB_NAME == *"gnu"* ]];   then COMPILER=GNU; fi
if [[ $JOB_NAME == *"intel"* ]]; then COMPILER=INTEL; fi
if [[ $JOB_NAME == *"cuda"* ]];  then COMPILER=CUDA; fi
if [[ $JOB_NAME == *"clang"* ]]; then COMPILER=CLANG; fi


# Set the optimization level
# Defaults to debug
if [[ $JOB_NAME == *"opt"* ]]; then export BUILD_TYPE=RELEASE; fi

# Set the node types default to serial
NODE_TYPE=SERIAL
if [[ $JOB_NAME == *"openmp"* ]]; then export USE_OPENMP=ON;   NODE_TYPE=OPENMP; fi
if [[ $JOB_NAME == *"cuda"* ]];   then export USE_CUDA=ON;     NODE_TYPE=CUDA;   fi
if [[ $JOB_NAME == *"pthread"* ]];then export USE_PTHREADS=ON; NODE_TYPE=THREAD; fi
if [[ $JOB_NAME == *"serial"* ]];then  NODE_TYPE=SERIAL; fi

if [[ $USE_OPENMP == "ON" && $USE_CUDA == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $USE_OPENMP == "ON" && $USE_PTHREADS == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $USE_OPENMP == "ON" && $NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $USE_CUDA == "ON" && $USE_PTHREADS == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $USE_CUDA == "ON" && $NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $USE_PTHREADS == "ON" && $NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
