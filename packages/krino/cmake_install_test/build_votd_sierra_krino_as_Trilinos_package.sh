#!/bin/sh

function execute() {
  stdbuf -o0 -e0 echo "% $@" ;
  if [ $# -gt 0 ] ; then
    if [ $1 == "module" ] ; then
      module "${@:2}" ;
    elif [ $1 == "source" ] ; then
      source "${@:2}" ;
    else
      eval "$@" ;
    fi
  fi
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

function build_yaml()
{
    productName=yaml
    buildType=$1

    cd_to_new_dir ${output_dir}/${productName}
    execute tar -xzf ${sierra_proj}/TPLs_src/spack/spack_tpls/yaml-cpp/*.tar.gz
    cd_to_new_dir ${productName}_build

    execute cmake -DCMAKE_BUILD_TYPE=${buildType} -DYAML_CPP_BUILD_TESTS=false -DCMAKE_INSTALL_PREFIX=../${productName}_install ../yaml-cpp
    make_and_install $productName
}

function setup_trilinos_with_krino()
{
    productName=trilinos
    
    cd_to_new_dir ${output_dir}/${productName}
    
    if [ ! -d Trilinos ] ; then
      execute git clone -b develop https://github.com/trilinos/Trilinos.git Trilinos
    else
      execute cd Trilinos
      execute git checkout develop
      execute git clean -fd
      execute git reset --hard
      execute git pull
      execute cd ..
    fi
    
    if [ -d Trilinos/packages/krino ] ; then
      execute rm -rf Trilinos/packages/krino;
      execute rm -rf Trilinos/packages/stk;
    fi
    if [ ! -L Trilinos/packages/krino ] ; then
      execute ln -s ${sierra_proj}/krino Trilinos/packages
      execute ln -s ${sierra_proj}/stk Trilinos/packages
    fi   
}

function build_trilinos_with_krino()
{
    productName=trilinos
    buildType=$1

    execute cd ${output_dir}/${productName}
    
    execute rm -rf ${productName}_${buildType}_install
    cd_to_new_dir ${productName}_${buildType}_build
    execute "TRILINOS_INSTALL_DIR=../${productName}_${buildType}_install CMAKE_BUILD_TYPE=${buildType} $sierra_proj/krino/cmake_install_test/run_cmake_krino"
    make_and_install $productName  
}

function setup_environment()
{
    execute source /etc/profile.d/modules.sh
    execute module purge
    
    execute source $sierra_proj/krino/cmake_install_test/load_gcc_modules

    execute module list
    execute env
}

function runTests()
{
    buildType=$1
    execute cd ${output_dir}/${productName}/${productName}_${buildType}_build
    execute ctest -j 16 || exit_with_message "$2 tests failed"
}

sierra_proj=${SIERRA_PROJ:-${PWD}}
output_dir=${OUTPUT_DIR:-${PWD}/../krino-cmake-testing}

cuda_on_or_off=${CUDA:-OFF}
date_suffix=`date +%F_%H-%M-%S`

if [ ! -d ${output_dir} ] ; then
  execute mkdir ${output_dir}
fi

setup_trilinos_with_krino
setup_environment

execute cd ${output_dir}

build_trilinos_with_krino release
build_yaml release
build_trilinos_with_krino release
runTests release
build_trilinos_with_krino debug
runTests debug
