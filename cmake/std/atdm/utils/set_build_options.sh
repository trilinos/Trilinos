# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Assert that JOB_NAME is set!
if [ -z "$JOB_NAME" ] ; then
  echo "Error, must set JOB_NAME in env!"
  return
fi

echo "Setting default compiler and build options for JOB_NAME='${JOB_NAME}'"

# Set the defaults
export ATDM_CONFIG_COMPILER=DEFAULT
export ATDM_CONFIG_BUILD_TYPE=DEBUG
export ATDM_CONFIG_USE_OPENMP=OFF
export ATDM_CONFIG_USE_CUDA=OFF
export ATDM_CONFIG_USE_PTHREADS=OFF

# Set the compiler
if [[ $JOB_NAME == "default" ]]; then
  ATDM_CONFIG_COMPILER=DEFAULT
elif [[ $JOB_NAME == *"gnu"* ]]; then
  ATDM_CONFIG_COMPILER=GNU
elif [[ $JOB_NAME == *"intel"* ]]; then
 ATDM_CONFIG_COMPILER=INTEL
elif [[ $JOB_NAME == *"cuda-8.0"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA-8.0
elif [[ $JOB_NAME == *"cuda-9.0"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA-9.0
elif [[ $JOB_NAME == *"cuda"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA
elif [[ $JOB_NAME == *"clang"* ]]; then
  ATDM_CONFIG_COMPILER=CLANG
else
  echo "***"
  echo "*** ERROR: A compiler was not specified in '$JOB_NAME'!"
  echo "***"
fi

# Set the optimization level
# Defaults to debug
if [[ $JOB_NAME == *"opt"* ]]; then export ATDM_CONFIG_BUILD_TYPE=RELEASE; fi

# Set the node types default to serial
ATDM_CONFIG_NODE_TYPE=SERIAL
if [[ $JOB_NAME == *"openmp"* ]]; then export ATDM_CONFIG_USE_OPENMP=ON;   ATDM_CONFIG_NODE_TYPE=OPENMP; fi
if [[ $JOB_NAME == *"cuda"* ]];   then export ATDM_CONFIG_USE_CUDA=ON;     ATDM_CONFIG_NODE_TYPE=CUDA;   fi
if [[ $JOB_NAME == *"pthread"* ]];then export ATDM_CONFIG_USE_PTHREADS=ON; ATDM_CONFIG_NODE_TYPE=THREAD; fi
if [[ $JOB_NAME == *"serial"* ]];then  ATDM_CONFIG_NODE_TYPE=SERIAL; fi

# Assert a consistent set of options
if [[ $ATDM_CONFIG_USE_OPENMP == "ON" && $ATDM_CONFIG_USE_CUDA == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_OPENMP == "ON" && $ATDM_CONFIG_USE_PTHREADS == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_OPENMP == "ON" && $ATDM_CONFIG_NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_CUDA == "ON" && $ATDM_CONFIG_USE_PTHREADS == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_CUDA == "ON" && $ATDM_CONFIG_NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_PTHREADS == "ON" && $ATDM_CONFIG_NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
