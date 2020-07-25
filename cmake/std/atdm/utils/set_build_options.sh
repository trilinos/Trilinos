#
# Script to parse a <build-name> string and set some env vars
#
# Usage:
#
#   source set_build_options.sh
#
# where ATDM_CONFIG_BULID_NAME must be set in the env prior to sourcing this
# scirpt.
#
# On completion, this will set the env vars:
#
#   ATDM_CONFIG_COMPILER
#   ATDM_CONFIG_USE_MPI
#   ATDM_CONFIG_KOKKOS_ARCH
#   ATDM_CONFIG_BUILD_TYPE
#   ATDM_CONFIG_NODE_TYPE
#   ATDM_CONFIG_USE_OPENMP
#   ATDM_CONFIG_USE_CUDA
#   ATDM_CONFIG_USE_PTHREADS
#   ATDM_CONFIG_CUDA_RDC
#   ATDM_CONFIG_FPIC
#   ATDM_CONFIG_COMPLEX
#   ATDM_CONFIG_SHARED_LIBS
#   ATDM_CONFIG_PT_PACKAGES
#
# or will error out.
#
# NOTE: Before sourcing this script consider sourcing
#
#   source unset_atdm_config_vars_build_options.sh
#
# or even:
#
#   source unset_atdm_config_vars_all.sh
#
# to wipe out old values of these vars, just to be safe.
#

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi
unset called

if [ -z "$ATDM_CONFIG_BUILD_NAME" ] ; then
  echo "Error, must set ATDM_CONFIG_BUILD_NAME in env!"
  return
fi

ATDM_UTILS_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
export ATDM_CONFIG_SCRIPT_DIR=`readlink -f ${ATDM_UTILS_SCRIPT_DIR}/..`

. ${ATDM_CONFIG_SCRIPT_DIR}/utils/define_atdm_match_keyword.sh

echo "Setting compiler and build options for build-name '${ATDM_CONFIG_BUILD_NAME}'"

# Set the defaults
export ATDM_CONFIG_COMPILER=DEFAULT
export ATDM_CONFIG_USE_MPI=ON
export ATDM_CONFIG_KOKKOS_ARCH=DEFAULT
export ATDM_CONFIG_BUILD_TYPE=DEBUG
export ATDM_CONFIG_SHARED_LIBS=OFF
export ATDM_CONFIG_NODE_TYPE=SERIAL
export ATDM_CONFIG_USE_OPENMP=OFF
export ATDM_CONFIG_USE_CUDA=OFF
export ATDM_CONFIG_USE_PTHREADS=OFF
export ATDM_CONFIG_CUDA_RDC=OFF
export ATDM_CONFIG_FPIC=OFF
export ATDM_CONFIG_COMPLEX=OFF
export ATDM_CONFIG_SHARED_LIBS=OFF
export ATDM_CONFIG_PT_PACKAGES=OFF

# Process system custom build logic
export ATDM_CONFIG_CUSTOM_COMPILER_SET=0
if [ -e ${ATDM_CONFIG_SYSTEM_DIR}/custom_builds.sh ]; then
  source ${ATDM_CONFIG_SYSTEM_DIR}/custom_builds.sh
fi

# NOTE: Currently only the specialization of ATDM_CONFIG_COMPILER from
# custom_builds.sh is supported below.  ToDo: Add support for customizing
# other things as needed.

# Set ATDM_CONFIG_COMPILER
if [[ "${ATDM_CONFIG_COMPILER}" != "DEFAULT" ]] ; then
  # Custom compile already set
  export ATDM_CONFIG_CUSTOM_COMPILER_SET=1
elif atdm_match_buildname_keyword default; then
  export ATDM_CONFIG_COMPILER=DEFAULT
elif atdm_match_buildname_keyword cuda-8.0; then
  export ATDM_CONFIG_COMPILER=CUDA-8.0
elif atdm_match_buildname_keyword cuda-9.0; then
  export ATDM_CONFIG_COMPILER=CUDA-9.0
elif atdm_match_any_buildname_keyword cuda-9.2-gnu-7.2.0 cuda-9.2_gnu-7.2.0 \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-9.2_GNU-7.2.0
elif atdm_match_buildname_keyword cuda-9.2; then
  export ATDM_CONFIG_COMPILER=CUDA-9.2
elif atdm_match_any_buildname_keyword cuda-10.0-gnu-7.4.0 cuda-10.0_gnu-7.4.0 \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.0_GNU-7.4.0
elif atdm_match_any_buildname_keyword cuda-10.1-gnu-7.2.0 cuda-10.1_gnu-7.2.0 \
  ; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1_GNU-7.2.0
elif atdm_match_buildname_keyword cuda-10.1; then
  export ATDM_CONFIG_COMPILER=CUDA-10.1
elif atdm_match_buildname_keyword cuda-10; then
  export ATDM_CONFIG_COMPILER=CUDA-10
elif atdm_match_buildname_keyword cuda; then
  export ATDM_CONFIG_COMPILER=CUDA
elif atdm_match_buildname_keyword gnu-4.8.4; then
  export ATDM_CONFIG_COMPILER=GNU-4.8.4
elif atdm_match_buildname_keyword gnu-4.9.3; then
  export ATDM_CONFIG_COMPILER=GNU-4.9.3
elif atdm_match_buildname_keyword gnu-6.1.0; then
  export ATDM_CONFIG_COMPILER=GNU-6.1.0
elif atdm_match_buildname_keyword gnu-7.2.0; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0
elif atdm_match_buildname_keyword gnu-7.4.0; then
  export ATDM_CONFIG_COMPILER=GNU-7.4.0
elif atdm_match_buildname_keyword gnu; then
  export ATDM_CONFIG_COMPILER=GNU
elif atdm_match_buildname_keyword intel-17.0.1; then
 export ATDM_CONFIG_COMPILER=INTEL-17.0.1
elif atdm_match_buildname_keyword intel-17; then
 export ATDM_CONFIG_COMPILER=INTEL-17.0.1
elif atdm_match_buildname_keyword intel-18.0.2; then
 export ATDM_CONFIG_COMPILER=INTEL-18.0.2
elif atdm_match_buildname_keyword intel-18.0.5; then
 export ATDM_CONFIG_COMPILER=INTEL-18.0.5
elif atdm_match_buildname_keyword intel-18; then
 export ATDM_CONFIG_COMPILER=INTEL-18.0.5
elif atdm_match_buildname_keyword intel; then
 export ATDM_CONFIG_COMPILER=INTEL
elif atdm_match_buildname_keyword clang-3.9.0; then
  export ATDM_CONFIG_COMPILER=CLANG-3.9.0
elif atdm_match_buildname_keyword clang-5.0.1; then
  export ATDM_CONFIG_COMPILER=CLANG-5.0.1
elif atdm_match_buildname_keyword clang-7.0.1; then
  export ATDM_CONFIG_COMPILER=CLANG-7.0.1
elif atdm_match_buildname_keyword clang; then
  export ATDM_CONFIG_COMPILER=CLANG
else
  echo
  echo "***"
  echo "*** ERROR: A compiler was not specified in '$ATDM_CONFIG_BUILD_NAME'!"
  echo "***"
fi
# NOTE: Above, the 'cuda' keywords need to be parsed first since they could
# have the compiler keywords embedded in them.  For example we need to match
# 'cuda-10.0-gnu-7.4.0' before we match 'gnu-7.4.0'.

# Use MPI or not
if   atdm_match_buildname_keyword no-mpi; then
  export ATDM_CONFIG_USE_MPI=OFF
elif atdm_match_buildname_keyword mpi; then
  export ATDM_CONFIG_USE_MPI=ON
fi

# Set ATDM_CONFIG_KOKKOS_ARCH
export ATDM_CONFIG_KOKKOS_ARCH=DEFAULT
. ${ATDM_UTILS_SCRIPT_DIR}/kokkos_arch_array.sh
for kokkos_arch in ${kokkos_arch_array[@]} ; do
  if atdm_match_buildname_keyword ${kokkos_arch}; then
    export ATDM_CONFIG_KOKKOS_ARCH=${kokkos_arch}
    break
  fi
done
if   [[ $ATDM_CONFIG_KOKKOS_ARCH == "DEFAULT" ]] \
  && [[ $ATDM_CONFIG_VERBOSE == "1" ]] \
   ; then
  echo "No KOKKOS_ARCH specified so using system default"
fi

# Set ATDM_CONFIG_BUILD_TYPE
if atdm_match_any_buildname_keyword \
    release-debug release_debug opt-dbg opt_dbg \
  ; then
  export ATDM_CONFIG_BUILD_TYPE=RELEASE-DEBUG
elif atdm_match_any_buildname_keyword \
  release opt \
  ; then
  export ATDM_CONFIG_BUILD_TYPE=RELEASE
elif atdm_match_any_buildname_keyword \
  debug dbg \
  ; then
  export ATDM_CONFIG_BUILD_TYPE=DEBUG
fi

# Set the node type vars
if   atdm_match_buildname_keyword cuda; then
  export ATDM_CONFIG_USE_CUDA=ON
  export ATDM_CONFIG_NODE_TYPE=CUDA
elif atdm_match_buildname_keyword openmp; then
  export ATDM_CONFIG_USE_OPENMP=ON
  export ATDM_CONFIG_NODE_TYPE=OPENMP
elif atdm_match_buildname_keyword pthread; then
  echo
  echo "***"
  echo "*** ERROR: The Kokkos Pthreads backend is no longer supported (see TRIL-272)!"
  echo "*** Please use a different backend like 'serial', 'openmp', or 'cuda'."
  echo "***"
  echo
  return
elif atdm_match_buildname_keyword serial; then
  export ATDM_CONFIG_NODE_TYPE=SERIAL
fi

# Use CUDA RDC or not
if   atdm_match_buildname_keyword no-rdc; then
  export ATDM_CONFIG_CUDA_RDC=OFF
elif atdm_match_buildname_keyword rdc; then
  export ATDM_CONFIG_CUDA_RDC=ON
fi

# Use -fPIC or not
if atdm_match_buildname_keyword fpic; then
  export ATDM_CONFIG_FPIC=ON
fi

# Enable complex (double) data-types or not
if atdm_match_buildname_keyword no-complex; then
  export ATDM_CONFIG_COMPLEX=OFF
elif atdm_match_buildname_keyword complex; then
  export ATDM_CONFIG_COMPLEX=ON
fi

# Set 'static' or 'shared'
if atdm_match_buildname_keyword shared; then
  export ATDM_CONFIG_SHARED_LIBS=ON
elif atdm_match_buildname_keyword static; then
  export ATDM_CONFIG_SHARED_LIBS=OFF
fi

# Allow enable of all Primary Tested (pt) packages are not
if atdm_match_buildname_keyword pt; then
  export ATDM_CONFIG_PT_PACKAGES=ON
fi

export ATDM_CONFIG_FINISHED_SET_BUILD_OPTIONS=1
