################################################################################
#
# Set up env on ats2 for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

# ats2 jobs all use the same environmnet changes to the
# sourced script below will impact jobs on both of those
# machines. please be mindful of this when making changes

# TODO: Check with Ross that pwr9 should be default
# TODO: Should pwr9 be in the build name?
if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  export ATDM_CONFIG_KOKKOS_ARCH=Power9
#TODO:  export ATDM_CONFIG_KOKKOS_ARCH=Volta70?
#TODO:  export ATDM_CONFIG_KOKKOS_ARCH=Volta72?
fi

if [ "$ATDM_CONFIG_KOKKOS_ARCH" != "Power9" ] ; then
  echo "***"
  echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option on this system."
  echo "*** '$ATDM_CONFIG_KOKKOS_ARCH' appears in $ATDM_CONFIG_BUILD_NAME which then sets the KOKKOS_ARCH"
  echo "*** on ats2 'Power9, Volta70?, Volta72?' are the only valid KOKKOS_ARCH. If no KOKKOS_ARCH is specified then"
  echo "*** 'Power9' will be used by default"
  echo "***"
  return
fi

export ATDM_CONFIG_SPARC_TPL_BASE=/projects/sparc/tpls/ats2-pwr9

export ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT=4:00:00

export ATDM_CONFIG_SLURM_DEFAULT_ACCOUNT=fy150090

source $ATDM_SCRIPT_DIR/common/toss3/environment_new.sh

export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=/projects/atdm_devops/trilinos_installs/
