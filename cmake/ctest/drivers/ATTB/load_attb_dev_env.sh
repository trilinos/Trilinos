if [[ `hostname` =~ hansen.* ]] ; then
  echo "Is hansen!"
  export ATTB_DEVPACK_MODULE=devpack/openmpi/1.10.0/gcc/4.8.4/cuda/7.5.18
elif [[ `hostname` =~ shiller.* ]] ; then
  echo "Is shiller!"
  export ATTB_DEVPACK_MODULE=devpack/openmpi/1.10.0/gcc/4.8.4/cuda/7.5.18
elif [[ `hostname` =~ shepard.* ]] ; then
  echo "Is shepard!"
  export ATTB_DEVPACK_MODULE=devpack/openmpi/1.10.0/intel/16.1.056/cuda/none
else
  echo "On some unknown machine, aborting!" 
fi

echo
echo "Setting up default dev env $ATTB_DEVPACK_MODULE"
module load $ATTB_DEVPACK_MODULE

echo
echo "Loading more recent git"
module load git/20150310
which git
git --version

echo
echo "Loading CMake 3.4.3"
module load cmake/3.4.3
which cmake
cmake --version
