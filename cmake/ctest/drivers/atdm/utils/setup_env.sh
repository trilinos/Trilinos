source ${WORKSPACE}/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME
echo
module list
echo
echo "cmake in path:"
which cmake
echo
echo "ninja in path:"
which ninja
echo
echo "ATDM config env vars:"
set | grep ATDM_CONFIG_
echo
echo "PATH=$PATH"

if [ "${Trilinos_REPOSITORY_LOCATION}" == "" ] ; then
  export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
fi
