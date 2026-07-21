#!/bin/bash

exe() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed.";
    exit 1;
  fi
}

output_dir=${OUTPUT_DIR:-/fgs/$USER/stk-cmake-testing}

if [[ -d "${output_dir}/spack" ]] ; then
  exe source ${output_dir}/spack/share/spack/setup-env.sh
  exe spack env activate stk_gcc12_spack_env
  exe module load aue/gcc/12.3.0
  exe g++ --version
else
  echo "Failed to find spack-dir ${output_dir}/spack. Run stk_spack_create_env.sh"
  return 1;
fi

return 0;
