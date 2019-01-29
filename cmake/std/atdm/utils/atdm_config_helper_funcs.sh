#
# This function aids in the building of an explicit list of libraries.
#
# Usage:
#
#   atdm_config_add_libs_to_var <var_name> <lib_dir> <ext> <lib0> <lib1> ...
#
# If the var <var_name> is already set, then it will be appended to.  On
# completion of the function, <var_name> is exported.  The result is a
# semi-colon seprated list of full library paths.  CMake likes this type of
# input for libraries and it makes them easier to verify that they exist.
#
# For example, to build up a set of libs for MLK, one would use something like:
# 
#   atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${CBLAS_ROOT}/mkl/lib/intel64 .so \
#      mkl_intel_lp64 mkl_intel_thread mkl_core
#
#   atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${CBLAS_ROOT}/lib/intel64 .so \
#      iomp5
#


function atdm_config_add_libs_to_var {
  
  # Formal arguments
  export_env_var_name=$1 ; shift
  lib_dir=$1 ; shift
  lib_ext=$1 ; shift
  libs_to_add=$@

  #echo
  #echo "export_env_var_name = '${export_env_var_name}'"
  #echo "lib_dir = '${lib_dir}'"
  #echo "lib_ext = '${lib_ext}'"
  #echo "libs_to_add = '${libs_to_add}'"

  full_libs_paths="${!export_env_var_name}"
  for lib_to_add in ${libs_to_add} ; do
    full_lib_path="${lib_dir}/lib${lib_to_add}${lib_ext}"
    if [[ "${full_libs_paths}" == "" ]] ; then
      full_libs_paths="${full_lib_path}"
    else
      full_libs_paths="${full_libs_paths};${full_lib_path}"
    fi
  done

  export ${export_env_var_name}=${full_libs_paths}
  #echo "${export_env_var_name} = '${!export_env_var_name}'"
  #echo

  unset export_env_var_name

}
