#!/bin/bash
################################################################################
# @Brief: On the specified arch, build and run xpy.
#
# Author: Evan Harvey <eharvey@sandia.gov>
################################################################################

function envprint() {
  for x in $@; do
      echo $x:\$$x | envsubst
  done
}

function printhelp() {
  echo "--Usage--"
  echo "$0 HOST_ARCH <ACCELERATOR_ARCH>"
  echo "  HOST_ARCH:        POWER9, A64FX, SKX"
  echo "  ACCELERATOR_ARCH: VOLTA70, AMPERE80"
  echo ""
  echo "Invocations used to collect us-rse-escience-2022 results:"
  echo "env KOKKOS_SRC_DIR=$HOME/KOKKOS.base/kokkos KOKKOSKERNELS_SRC_DIR=$HOME/KOKKOS.base/kokkos-kernels/ KOKKOSKERNELS_SHA=TODO-HEAD-SHA ./us-rse-escience-2022-reproducer.sh POWER9 VOLTA70"
  echo "env KOKKOS_SRC_DIR=$HOME/KOKKOS.base/kokkos KOKKOSKERNELS_SRC_DIR=$HOME/KOKKOS.base/kokkos-kernels/ KOKKOSKERNELS_SHA=TODO-HEAD-SHA ./us-rse-escience-2022-reproducer.sh AMPERE80"
}

function earlyexit() {
   rm -rf $benchmark_dir
   exit $1
}

function beval() {
  local ret=0
  echo "---------------------------------------------------------------------------------------------------------------"
  echo "START: \"$@\""
  if [ $dry_run == "off" ]; then
    eval $@
    ret=$PIPESTATUS
  fi
  if [ $ret -ne 0 ]; then
      echo "ERROR: \"$@\""
      earlyexit 1
  fi
  echo "END  : \"$@\""
  echo "---------------------------------------------------------------------------------------------------------------"
}

# Handle input args
export KOKKOS_SRC_DIR=${KOKKOS_SRC_DIR:-"$HOME/KOKKOS.base/kokkos"}
export KOKKOS_SRC_DIR=$(realpath $KOKKOS_SRC_DIR)
export KOKKOS_SHA=${KOKKOS_SHA:-"tags/3.6.00"}
export KOKKOSKERNELS_SRC_DIR=${KOKKOSKERNELS_SRC_DIR:-"$HOME/KOKKOS.base/kokkos-kernels"}
export KOKKOSKERNELS_SRC_DIR=$(realpath $KOKKOSKERNELS_SRC_DIR)
envprint KOKKOS_SRC_DIR KOKKOS_SHA KOKKOSKERNELS_SRC_DIR KOKKOSKERNELS_SHA

dry_run="off"
arch_names="$1 $2"
echo "HOST_ARCH=\"$1\", ACCELERATOR_ARCH=\"$2\""

# Create benchmark directory
benchmark_dir=$PWD/$0_$(date +"%Y-%m-%d_%H.%M.%S")
beval mkdir -p $benchmark_dir/kokkos-{build,install}
beval mkdir -p $benchmark_dir/kokkos-kernels-{build,install}
export KOKKOS_BUILD_DIR=$(realpath $benchmark_dir/kokkos-build)
export KOKKOS_INSTALL_DIR=$(realpath $benchmark_dir/kokkos-install)
export KOKKOSKERNELS_BUILD_DIR=$(realpath $benchmark_dir/kokkos-kernels-build)
export KOKKOSKERNELS_INSTALL_DIR=$(realpath $benchmark_dir/kokkos-kernels-install)
envprint KOKKOS_INSTALL_DIR KOKKOS_BUILD_DIR KOKKOSKERNELS_BUILD_DIR KOKKOSKERNELS_INSTALL_DIR

# Setup arch specific cmake configurations and job submission commands
if [[ "$arch_names" == " " ]]; then
    printhelp; earlyexit 1
elif [ "$arch_names" == "POWER9 VOLTA70" ]; then
  module purge
  module load cuda/11.2.0 gcc/8.3.1 cmake/3.18.0
  kokkos_config_cmd="cd $KOKKOS_BUILD_DIR && $KOKKOS_SRC_DIR/generate_makefile.bash --cxxflags='-O3' \
                     --arch=Power9,Volta70 --with-cuda=$CUDA_PATH --compiler=$KOKKOS_SRC_DIR/bin/nvcc_wrapper \
                     --kokkos-path=$KOKKOS_SRC_DIR --prefix=$KOKKOS_INSTALL_DIR 2>&1 | tee kokkos_config_cmd.out"
  kokkos_config_defaults_cmd="cd $KOKKOS_BUILD_DIR && cmake -DKokkos_ENABLE_TESTS:BOOL=OFF $KOKKOS_SRC_DIR 2>&1 \
                              | tee -a kokkos_config_cmd.out"

  kokkoskernels_config_cmd="cd $KOKKOSKERNELS_BUILD_DIR && $KOKKOSKERNELS_SRC_DIR/cm_generate_makefile.bash \
                            --arch=Power9,Volta70 --with-cuda=$CUDA_PATH -- --compiler=$KOKKOS_INSTALL_DIR/bin/nvcc_wrapper \
                            --cxxflags='-O3' --disable-tests --enable-examples \
                            --kokkos-path=$KOKKOS_SRC_DIR --kokkoskernels-path=$KOKKOSKERNELS_SRC_DIR \
                            --kokkos-prefix=$KOKKOS_INSTALL_DIR --prefix=$KOKKOSKERNELS_INSTALL_DIR 2>&1 | \
                            tee kokkoskernels_config_cmd.out"
  kokkoskernels_config_defaults_cmd="cd $KOKKOSKERNELS_BUILD_DIR && cmake -DKokkosKernels_INST_LAYOUTLEFT:BOOL=OFF \
                                   -DKokkosKernels_INST_LAYOUTRIGHT:BOOL=ON -DKokkosKernels_INST_DOUBLE:BOOL=OFF \
                                   $KOKKOSKERNELS_SRC_DIR 2>&1 | tee -a kokkoskernels_config_cmd.out"

  kokkos_build_cmd="bsub -q normal -W 2:00 -Is $KOKKOS_BUILD_DIR/build.sh"
  kokkoskernels_build_cmd="bsub -q normal -W 2:00 -Is $KOKKOSKERNELS_BUILD_DIR/build.sh"
  benchmark_cmd="bsub -q normal -W 2:00 -Is $KOKKOSKERNELS_BUILD_DIR/bench.sh"
elif [ "$arch_names" == "SNB VOLTA70" ]; then
  module purge
  module load sems-archive-env sems-env sems-gcc/8.3.0 sems-cmake/3.19.1 cuda/11.2 sems-archive-git/2.10.1
  kokkos_config_cmd="cd $KOKKOS_BUILD_DIR && $KOKKOS_SRC_DIR/generate_makefile.bash --cxxflags='-O3' \
                     --arch=SNB,Volta70 --with-cuda=$CUDA_PATH --compiler=$KOKKOS_SRC_DIR/bin/nvcc_wrapper \
                     --kokkos-path=$KOKKOS_SRC_DIR --prefix=$KOKKOS_INSTALL_DIR 2>&1 | tee kokkos_config_cmd.out"
  kokkos_config_defaults_cmd="cd $KOKKOS_BUILD_DIR && cmake -DKokkos_ENABLE_TESTS:BOOL=OFF $KOKKOS_SRC_DIR 2>&1 \
                              | tee -a kokkos_config_cmd.out"

  kokkoskernels_config_cmd="cd $KOKKOSKERNELS_BUILD_DIR && $KOKKOSKERNELS_SRC_DIR/cm_generate_makefile.bash \
                            --arch=SNB,Volta70 --with-cuda=$CUDA_PATH --compiler=$KOKKOS_INSTALL_DIR/bin/nvcc_wrapper \
                            --cxxflags='-O3' --disable-tests --enable-examples \
                            --kokkos-path=$KOKKOS_SRC_DIR --kokkoskernels-path=$KOKKOSKERNELS_SRC_DIR \
                            --kokkos-prefix=$KOKKOS_INSTALL_DIR --prefix=$KOKKOSKERNELS_INSTALL_DIR 2>&1 | \
                            tee kokkoskernels_config_cmd.out"
  kokkoskernels_config_defaults_cmd="cd $KOKKOSKERNELS_BUILD_DIR && cmake -DKokkosKernels_INST_LAYOUTLEFT:BOOL=OFF \
                                   -DKokkosKernels_INST_LAYOUTRIGHT:BOOL=ON -DKokkosKernels_INST_DOUBLE:BOOL=OFF \
                                   $KOKKOSKERNELS_SRC_DIR 2>&1 | tee -a kokkoskernels_config_cmd.out"

  kokkos_build_cmd="$KOKKOS_BUILD_DIR/build.sh"
  kokkoskernels_build_cmd="$KOKKOSKERNELS_BUILD_DIR/build.sh"
  benchmark_cmd="$KOKKOSKERNELS_BUILD_DIR/bench.sh"
elif [ "$arch_names" == "AMPERE80" ]; then
  module purge
  module load cudatoolkit/11.2 cmake/3.22.0

  kokkos_config_cmd="cd $KOKKOS_BUILD_DIR && $KOKKOS_SRC_DIR/generate_makefile.bash --cxxflags='-O3' \
                    --arch=Ampere80 --with-cuda=$CUDA_HOME --compiler=$KOKKOS_SRC_DIR/bin/nvcc_wrapper \
                    --kokkos-path=$KOKKOS_SRC_DIR --prefix=$KOKKOS_INSTALL_DIR &> kokkos_config_cmd.out"

  kokkos_config_defaults_cmd="cd $KOKKOS_BUILD_DIR && cmake -DKokkos_ENABLE_TESTS:BOOL=OFF $KOKKOS_SRC_DIR &>  kokkos_config_cmd.out"
  kokkoskernels_config_cmd="cd $KOKKOSKERNELS_BUILD_DIR && $KOKKOSKERNELS_SRC_DIR/cm_generate_makefile.bash \
                           --arch=Ampere80 --with-cuda=$CUDA_HOME --compiler=$KOKKOS_INSTALL_DIR/bin/nvcc_wrapper \
                           --cxxflags='-O3' --disable-tests --enable-examples \
                           --kokkos-path=$KOKKOS_SRC_DIR --kokkoskernels-path=$KOKKOSKERNELS_SRC_DIR \
                           --kokkos-prefix=$KOKKOS_INSTALL_DIR --prefix=$KOKKOSKERNELS_INSTALL_DIR &> kokkoskernels_config_cmd.out"

  kokkoskernels_config_defaults_cmd="cd $KOKKOSKERNELS_BUILD_DIR && cmake -S $KOKKOSKERNELS_SRC_DIR -DKokkosKernels_INST_LAYOUTLEFT:BOOL=OFF \
                                  -DKokkosKernels_INST_LAYOUTRIGHT:BOOL=ON -DKokkosKernels_INST_DOUBLE:BOOL=OFF &> kokkoskernels_config_cmd.out"

  kokkos_build_cmd="$KOKKOS_BUILD_DIR/build.sh"
  kokkoskernels_build_cmd="$KOKKOSKERNELS_BUILD_DIR/build.sh"
  benchmark_cmd="$KOKKOSKERNELS_BUILD_DIR/bench.sh"
elif [ "$arch_names" == "A64FX " ]; then
  export OMP_PROC_BIND=close
  export OMP_PLACES=cores
  export OMP_NUM_THREADS=48
  module purge
  module load gcc/10.2.0 cmake/3.17.0
  kokkos_config_cmd="cd $KOKKOS_BUILD_DIR && $KOKKOS_SRC_DIR/generate_makefile.bash --cxxflags='-O3' \
                     --arch=A64FX \
                     --kokkos-path=$KOKKOS_SRC_DIR --prefix=$KOKKOS_INSTALL_DIR 2>&1 | tee kokkos_config_cmd.out"
  kokkos_config_defaults_cmd="cd $KOKKOS_BUILD_DIR && cmake -DKokkos_ENABLE_TESTS:BOOL=OFF $KOKKOS_SRC_DIR 2>&1 \
                              | tee -a kokkos_config_cmd.out"

  kokkoskernels_config_cmd="cd $KOKKOSKERNELS_BUILD_DIR && $KOKKOSKERNELS_SRC_DIR/cm_generate_makefile.bash \
                            --cxxflags='-msve-vector-bits=512 -Ofast' --arch=A64FX --with-openmp \
                            --disable-tests --enable-examples \
                            --kokkos-path=$KOKKOS_SRC_DIR --kokkoskernels-path=$KOKKOSKERNELS_SRC_DIR \
                            --kokkos-prefix=$KOKKOS_INSTALL_DIR --prefix=$KOKKOSKERNELS_INSTALL_DIR 2>&1 | \
                            tee kokkoskernels_config_cmd.out"
  kokkoskernels_config_defaults_cmd="cd $KOKKOSKERNELS_BUILD_DIR && cmake -DKokkosKernels_INST_LAYOUTLEFT:BOOL=OFF \
                                   -DKokkosKernels_INST_LAYOUTRIGHT:BOOL=ON -DKokkosKernels_INST_DOUBLE:BOOL=OFF \
                                   $KOKKOSKERNELS_SRC_DIR 2>&1 | tee -a kokkoskernels_config_cmd.out"

  kokkos_build_cmd="srun --time=2:00:00 -N1 $KOKKOS_BUILD_DIR/build.sh"
  kokkoskernels_build_cmd="srun --time=2:00:00 -N1 $KOKKOSKERNELS_BUILD_DIR/build.sh"
  benchmark_cmd="srun --time=2:00:00 -N1 $KOKKOSKERNELS_BUILD_DIR/bench.sh"
elif [ "$arch_names" == "SKX " ]; then
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export OMP_NUM_THREADS=96
    module purge
    module load gcc/7.2.0 cmake/3.19.3
    kokkos_config_cmd="cd $KOKKOS_BUILD_DIR && $KOKKOS_SRC_DIR/generate_makefile.bash --cxxflags='-O3' \
                       --arch=SKX \
                       --kokkos-path=$KOKKOS_SRC_DIR --prefix=$KOKKOS_INSTALL_DIR 2>&1 | tee kokkos_config_cmd.out"
    kokkos_config_defaults_cmd="cd $KOKKOS_BUILD_DIR && cmake -DKokkos_ENABLE_TESTS:BOOL=OFF $KOKKOS_SRC_DIR 2>&1 \
                                | tee -a kokkos_config_cmd.out"

    kokkoskernels_config_cmd="cd $KOKKOSKERNELS_BUILD_DIR && $KOKKOSKERNELS_SRC_DIR/cm_generate_makefile.bash \
                              --cxxflags='-O3' --arch=SKX --with-openmp --disable-tests --enable-examples \
                              --kokkos-path=$KOKKOS_SRC_DIR --kokkoskernels-path=$KOKKOSKERNELS_SRC_DIR \
                              --kokkos-prefix=$KOKKOS_INSTALL_DIR --prefix=$KOKKOSKERNELS_INSTALL_DIR 2>&1 | \
                              tee kokkoskernels_config_cmd.out"
    kokkoskernels_config_defaults_cmd="cd $KOKKOSKERNELS_BUILD_DIR && cmake -DKokkosKernels_INST_LAYOUTLEFT:BOOL=OFF \
                                     -DKokkosKernels_INST_LAYOUTRIGHT:BOOL=ON -DKokkosKernels_INST_DOUBLE:BOOL=OFF \
                                     $KOKKOSKERNELS_SRC_DIR 2>&1 | tee -a kokkoskernels_config_cmd.out"

    kokkos_build_cmd="srun --time=2:00:00 -N1 $KOKKOS_BUILD_DIR/build.sh"
    kokkoskernels_build_cmd="srun --time=2:00:00 -N1 $KOKKOSKERNELS_BUILD_DIR/build.sh"
    benchmark_cmd="srun --time=2:00:00 -N1 $KOKKOSKERNELS_BUILD_DIR/bench.sh"
    use_simd="--use_simd=1"
else
  echo "Invalid arch: $arch_names"
  printhelp; earlyexit 1
fi

# Write the arch agnostic kokkos build script
echo "#!/bin/bash" > $KOKKOS_BUILD_DIR/build.sh
echo "cd $KOKKOS_BUILD_DIR" >> $KOKKOS_BUILD_DIR/build.sh
echo "make -j40 install" >> $KOKKOS_BUILD_DIR/build.sh
chmod +x $KOKKOS_BUILD_DIR/build.sh

# Write the arch agnostic kokkos-kernels build script
echo "#!/bin/bash" > $KOKKOSKERNELS_BUILD_DIR/build.sh
echo "cd $KOKKOSKERNELS_BUILD_DIR/example/half" >> $KOKKOSKERNELS_BUILD_DIR/build.sh
echo "make -j40 xpy" >> $KOKKOSKERNELS_BUILD_DIR/build.sh
chmod +x $KOKKOSKERNELS_BUILD_DIR/build.sh

# Write the arch agnostic kokkos-kernels benchmark script
echo "#!/bin/bash" > $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "cd $benchmark_dir" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 10 0 &> xpy_relative_error-10.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 100 0 &> xpy_relative_error-100.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 1000 0 &> xpy_relative_error-1000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 10000 0 &> xpy_relative_error-10000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 100000 0 &> xpy_relative_error-100000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh

echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 50000 1 &> xpy_runtime_only-50000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 500000 1 &> xpy_runtime_only-500000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 5000000 1 &> xpy_runtime_only-5000000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 50000000 1 &> xpy_runtime_only-50000000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
echo "$KOKKOSKERNELS_BUILD_DIR/example/half/xpy 500000000 1 &> xpy_runtime_only-500000000.out" >> $KOKKOSKERNELS_BUILD_DIR/bench.sh
chmod +x $KOKKOSKERNELS_BUILD_DIR/bench.sh

# Check out the correct SHAs
beval "cd $KOKKOS_SRC_DIR && git checkout $KOKKOS_SHA"
beval "cd $KOKKOSKERNELS_SRC_DIR && git checkout $KOKKOSKERNELS_SHA"

# Build Kokkos
beval $kokkos_config_cmd
beval $kokkos_config_defaults_cmd
beval $kokkos_build_cmd

# Wait for the file system on the head node to catch up
while [[ "$arch_names" == "POWER9 VOLTA70" && ! -e $KOKKOS_INSTALL_DIR/bin/nvcc_wrapper ]]; do
  sleep 3s
done

# Build KokkosKernels
beval $kokkoskernels_config_cmd
beval $kokkoskernels_config_defaults_cmd
beval $kokkoskernels_build_cmd

# Run the benchmark
beval $benchmark_cmd
