# This function aids in the building of an explicit list of libraries.
#
# Usage:
#
#   atdm_config_add_libs_to_var <var_name> <lib_dir> <ext> <lib0> <lib1> ...
#
# If the var <var_name> is already set, then it will be appended to.  On
# completion of the function, <var_name> is exported.  The result is a
# semi-colon separated list of full library paths.  CMake likes this type of
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


# Get the standard name of the sparc-dev module for a standard-named
# ATDM_CONFIG_COMPILER var.
#
function get_sparc_dev_module_name() {
  atdm_config_compiler=$1
  sparc_module_name=sparc-dev/$(echo "$ATDM_CONFIG_COMPILER" | tr '[:upper:]' '[:lower:]' | sed 's|/gnu-|/gcc-|g')
  #echo "sparc_module_name = '${sparc_module_name}'"
  sparc_module_name=$(echo "$sparc_module_name" | sed 's|gnu-|gcc-|g')
  #echo "sparc_module_name = '${sparc_module_name}'"
  echo "${sparc_module_name}"
}


# Remove the substrings from the environment variable.
#
# Usage:
#
#   atdm_remove_substrings_from_env_var <env_var> <delim> <sstr1> <sstr2> ...
#
# @param env_var:  The environment variable to modify.
# @param delim:    The delimiter used in the environment variable.
# @param sub_strs: The substrings to remove.
# @return void, the environment variable is exported to the new value.
#
function atdm_remove_substrings_from_env_var() {
  local env_var="$1"; shift
  local delim="$1"; shift
  local sub_strs="$@"
  #echo "${env_var}=${!env_var}"

  if [[ ! "$delim" = ":" &&
	      ! "$delim" = "," &&
        ! "$delim" = "-" &&
        ! "$delim" = "_" ]]; then
    printf "%s\n" "ERROR: $FUNCNAME: \"$delim\" is an invalid delimiter." 2>&1
    return
  fi

  local env_var_sub_strs=$(printf "%s" "${!env_var}" | sed "s/${delim}/ /g")

  for str in $env_var_sub_strs; do
    for subStr in $sub_strs; do
      if [ "$subStr" = "$str" ]; then
        env_var_sub_strs=("${env_var_sub_strs[@]/" $str"}")
        env_var_sub_strs=("${env_var_sub_strs[@]/"$str "}")
      fi
    done
  done

  export ${env_var}=$(printf "%s" "$env_var_sub_strs" | sed "s/ /$delim/g")
  #echo "${env_var}=${!env_var}"
}


# Remove individual dirs from the PATH environment variable
#
# Usage:
#
#   atdm_remove_dirs_from_path <dir1> <dir2> ...
#
# @param dirs: one or more space delimited individual directories
# @return void, the environment variable is exported to the new value.
#
function atdm_remove_dirs_from_path() {
  atdm_remove_substrings_from_env_var PATH ":" $@
}
