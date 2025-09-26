#!/bin/bash

exe() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed.";
    exit 1;
  fi
}

trilinos_dir=${1:?"Must pass directory to clone Trilinos into!"}
sierra_proj=${2:?"Must pass path to sierra!"}

if [ ! -d ${trilinos_dir} ] ; then
  exe git clone -b develop https://github.com/trilinos/Trilinos.git ${trilinos_dir}
else
  exe cd "${trilinos_dir}"
  exe git checkout develop
  exe git reset --hard origin/develop
  exe git pull
fi

if [ -d ${trilinos_dir}/packages/stk ] ; then
  exe rm -rf ${trilinos_dir}/packages/stk;
fi
if [ ! -L ${trilinos_dir}/packages/stk ] ; then
  exe ln -s ${sierra_proj}/stk ${trilinos_dir}/packages
fi

# we do not want to pay attention to the Trilinos configuration
#rm -rf ${trilinos_dir}/CTestConfig.cmake
