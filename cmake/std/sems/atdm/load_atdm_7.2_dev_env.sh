module purge

module load sems-env
module load atdm-env
module load sems-python/2.7.9
module load atdm-cmake/3.11.1
module load sems-git/2.10.1
module load atdm-gcc/7.2.0
module load atdm-openmpi/1.6.5/atdm
module load atdm-boost/1.63.0/atdm
module load atdm-zlib/1.2.8/atdm
module load atdm-hdf5/1.8.12/atdm
module load atdm-netcdf/4.4.1/atdm
module load atdm-parmetis/4.0.3/atdm
module load atdm-scotch/6.0.3/atdm
module load atdm-superlu/4.3/atdm

if [ "${TRILINOS_SEMS_DEV_ENV_VERBOSE}" == "1" ] ; then
  module list
fi

#
# D) Remember the loaded SEMS Dev Env
#

export TRILINOS_SEMS_DEV_ENV_LOADED=1
