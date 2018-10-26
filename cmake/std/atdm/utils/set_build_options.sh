# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Assert that ATDM_CONFIG_JOB_NAME is set!
if [ -z "$ATDM_CONFIG_JOB_NAME" ] ; then
  echo "Error, must set ATDM_CONFIG_JOB_NAME in env!"
  return
fi

echo "Setting default compiler and build options for ATDM_CONFIG_JOB_NAME='${ATDM_CONFIG_JOB_NAME}'"

# Set the defaults
export ATDM_CONFIG_COMPILER=DEFAULT
export ATDM_CONFIG_KOKKOS_ARCH=DEFAULT
export ATDM_CONFIG_BUILD_TYPE=DEBUG
export ATDM_CONFIG_USE_OPENMP=OFF
export ATDM_CONFIG_USE_CUDA=OFF
export ATDM_CONFIG_USE_PTHREADS=OFF

# Set the compiler
if [[ $ATDM_CONFIG_JOB_NAME == *"default" ]]; then
  ATDM_CONFIG_COMPILER=DEFAULT
elif [[ $ATDM_CONFIG_JOB_NAME == *"gnu-4.9.3"* ]]; then
  ATDM_CONFIG_COMPILER=GNU-4.9.3
elif [[ $ATDM_CONFIG_JOB_NAME == *"gnu-7.2.0"* ]]; then
  ATDM_CONFIG_COMPILER=GNU-7.2.0
elif [[ $ATDM_CONFIG_JOB_NAME == *"gnu"* ]]; then
  ATDM_CONFIG_COMPILER=GNU
elif [[ $ATDM_CONFIG_JOB_NAME == *"intel"* ]]; then
 ATDM_CONFIG_COMPILER=INTEL
elif [[ $ATDM_CONFIG_JOB_NAME == *"cuda-8.0"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA-8.0
elif [[ $ATDM_CONFIG_JOB_NAME == *"cuda-9.0"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA-9.0
elif [[ $ATDM_CONFIG_JOB_NAME == *"cuda-9.2"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA-9.2
elif [[ $ATDM_CONFIG_JOB_NAME == *"cuda"* ]]; then
  ATDM_CONFIG_COMPILER=CUDA
elif [[ $ATDM_CONFIG_JOB_NAME == *"clang"* ]]; then
  ATDM_CONFIG_COMPILER=CLANG
else
  echo
  echo "***"
  echo "*** ERROR: A compiler was not specified in '$ATDM_CONFIG_JOB_NAME'!"
  echo "***"
fi

# Set the KOKKOS_ARCH
if [[ $ATDM_CONFIG_JOB_NAME == *"-AMDAVX"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=AMDAVX
elif [[ $ATDM_CONFIG_JOB_NAME == *"-ARMv8-ThunderX"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=ARMv8-ThunderX
elif [[ $ATDM_CONFIG_JOB_NAME == *"-ARMv80"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=ARMv80
elif [[ $ATDM_CONFIG_JOB_NAME == *"-ARMv81"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=ARMv81
elif [[ $ATDM_CONFIG_JOB_NAME == *"-BDW"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=BDW
elif [[ $ATDM_CONFIG_JOB_NAME == *"-BGQ"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=BGQ
elif [[ $ATDM_CONFIG_JOB_NAME == *"-HSW"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=HSW
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Kepler30"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Kepler30
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Kepler32"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Kepler32
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Kepler35"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Kepler35
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Kepler37"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Kepler37
elif [[ $ATDM_CONFIG_JOB_NAME == *"-KNC"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=KNC
elif [[ $ATDM_CONFIG_JOB_NAME == *"-KNL"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=KNL
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Maxwell50"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Maxwell50
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Maxwell52"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Maxwell52
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Maxwell53"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Maxwell53
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Pascal60"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Pascal60
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Pascal61"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Pascal61
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Power7"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Power7
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Power8"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Power8
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Power9"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Power9
elif [[ $ATDM_CONFIG_JOB_NAME == *"-SKX"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=SKX
elif [[ $ATDM_CONFIG_JOB_NAME == *"-SNB"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=SNB
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Volta70"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Volta70
elif [[ $ATDM_CONFIG_JOB_NAME == *"-Volta72"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=Volta72
elif [[ $ATDM_CONFIG_JOB_NAME == *"-WSM"* ]]; then
  ATDM_CONFIG_KOKKOS_ARCH=WSM
else
  ATDM_CONFIG_KOKKOS_ARCH=DEFAULT
  echo "No KOKKOS_ARCH specified so using system default"
  # NOTE: <system_name>/environment.sh may set ATDM_CONFIG_KOKKOS_ARCH="" as
  # an allowed value!
fi

# Set the optimization level
# Defaults to debug
if [[ $ATDM_CONFIG_JOB_NAME == *"release-debug"* ]]; then
  export ATDM_CONFIG_BUILD_TYPE=RELEASE_DEBUG;
elif [[ $ATDM_CONFIG_JOB_NAME == *"opt"* ]]; then
  export ATDM_CONFIG_BUILD_TYPE=RELEASE;
elif [[ $ATDM_CONFIG_JOB_NAME == *"release"* ]]; then
  export ATDM_CONFIG_BUILD_TYPE=RELEASE;
elif [[ $ATDM_CONFIG_JOB_NAME == *"debug"* ]]; then
  export ATDM_CONFIG_BUILD_TYPE=DEBUG;
fi

# Set the node types default to serial
ATDM_CONFIG_NODE_TYPE=SERIAL
if [[ $ATDM_CONFIG_JOB_NAME == *"openmp"* ]]; then export ATDM_CONFIG_USE_OPENMP=ON;   ATDM_CONFIG_NODE_TYPE=OPENMP; fi
if [[ $ATDM_CONFIG_JOB_NAME == *"cuda"* ]];   then export ATDM_CONFIG_USE_CUDA=ON;     ATDM_CONFIG_NODE_TYPE=CUDA;   fi
if [[ $ATDM_CONFIG_JOB_NAME == *"pthread"* ]];then export ATDM_CONFIG_USE_PTHREADS=ON; ATDM_CONFIG_NODE_TYPE=THREAD; fi
if [[ $ATDM_CONFIG_JOB_NAME == *"serial"* ]];then  ATDM_CONFIG_NODE_TYPE=SERIAL; fi

# Assert a consistent set of options
if [[ $ATDM_CONFIG_USE_OPENMP == "ON" && $ATDM_CONFIG_USE_CUDA == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_OPENMP == "ON" && $ATDM_CONFIG_USE_PTHREADS == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_OPENMP == "ON" && $ATDM_CONFIG_NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_CUDA == "ON" && $ATDM_CONFIG_USE_PTHREADS == "ON" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_CUDA == "ON" && $ATDM_CONFIG_NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
if [[ $ATDM_CONFIG_USE_PTHREADS == "ON" && $ATDM_CONFIG_NODE_TYPE == "SERIAL" ]]; then echo "Can't set more than one backend"; return 1; fi
