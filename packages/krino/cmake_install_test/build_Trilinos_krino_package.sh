#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Please provide the path to the output directory as the only argument."
    exit 1
fi


function execute() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed.";
    exit 1;
  fi
}

function cd_to_new_dir()
{
    execute rm -rf $1
    execute mkdir $1
    execute cd $1
}

function exit_with_message()
{
    echo $1
    exit 1
}

function make_and_install()
{
    make -j 16 || exit_with_message "Failed building $1"
    make install || exit_with_message "Failed installing $1"
}

function pull_and_build_yaml()
{
    productName=yaml

    cd_to_new_dir ${output_dir}/${productName}
    
    execute git clone https://github.com/jbeder/yaml-cpp.git yaml-cpp
    cd yaml-cpp
    git checkout tags/yaml-cpp-0.6.3
    cd ..
    
    cd_to_new_dir ${productName}_build

    export CC=gcc
    export CXX=g++
    execute cmake -DCMAKE_BUILD_TYPE=${build_type^^} -DYAML_CPP_BUILD_TESTS=false -DCMAKE_INSTALL_PREFIX=../${productName}_install ../yaml-cpp
    make_and_install $productName
    unset CC
    unset CXX
}

function pull_Trilinos()
{
    productName=trilinos
    
    trilinos_dir=${output_dir}/${productName}
    if [ ! -d ${trilinos_dir} ] ; then
      execute mkdir ${trilinos_dir}
    fi
    execute cd ${trilinos_dir}
    
    if [ ! -d Trilinos ] ; then
      execute git clone -b develop https://github.com/trilinos/Trilinos.git Trilinos
    else
      execute cd Trilinos
      execute git checkout develop
      execute git reset --hard origin/develop
      execute git pull
      execute cd ..
    fi
}

function build_krino_Trilinos_package()
{
    productName=trilinos

    execute cd ${output_dir}/${productName}
    
    rm -rf ${productName}_install
    cd_to_new_dir ${productName}_build

    export TRILINOS_INSTALL_DIR=../${productName}_install
    
    execute ../Trilinos/packages/krino/cmake_install_test/run_cmake_krino
    make_and_install $productName    
}


function setup_environment()
{
    module purge
    
    # compiler environment
    module load cde/dev/cmake/3.19.2
    module load cde/dev/compiler/gcc/7.2.0
    module load cde/dev/gcc/7.2.0/openmpi/4.0.5
    
    #TPLs
    module load cde/dev/gcc/7.2.0/netlib-lapack/3.8.0
    module load cde/dev/gcc/7.2.0/hdf5/1.10.6
    module load cde/dev/gcc/7.2.0/netcdf-c/4.7.3
    module load cde/dev/gcc/7.2.0/parallel-netcdf/1.12.1
    module load cde/dev/gcc/7.2.0/metis/5.1.0
    module load cde/dev/gcc/7.2.0/parmetis/4.0.3
}

output_dir=$1

build_type=${CMAKE_BUILD_TYPE:-release}
date_suffix=`date +%F_%H-%M-%S`

if [ ! -d ${output_dir} ] ; then
  execute mkdir ${output_dir}
fi

pull_Trilinos
setup_environment

pull_and_build_yaml
build_krino_Trilinos_package
