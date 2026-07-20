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
reset_trilinos=${3:?"true"}

if [ ! -d ${trilinos_dir} ] ; then
  exe git clone -b develop https://github.com/trilinos/Trilinos.git ${trilinos_dir}
else
  if [[ ${reset_trilinos} == "true" ]] ; then
    exe cd "${trilinos_dir}"
    exe git checkout develop
    exe git reset --hard origin/develop
    exe git pull
  fi
fi

if [[ -L ${trilinos_dir}/packages/stk ]] ; then
  echo "${trilinos_dir}/packages/stk is sym-link..."
else
  exe rm -rf ${trilinos_dir}/packages/stk
  exe ln -s ${sierra_proj}/stk ${trilinos_dir}/packages
fi

