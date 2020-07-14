################################################################################
#
# Set up env on cts1 for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

# Chama and cts1 jobs all use the same environmnet changes to the
# sourced script below will impact jobs on both of those
# machines. please be mindful of this when making changes

if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  export ATDM_CONFIG_KOKKOS_ARCH=BDW
fi

if [ "$ATDM_CONFIG_KOKKOS_ARCH" != "BDW" ] ; then
  echo "***"
  echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option on this system."
  echo "*** '$ATDM_CONFIG_KOKKOS_ARCH' appears in $ATDM_CONFIG_BUILD_NAME which then sets the KOKKOS_ARCH"
  echo "*** on cts1 'BDW' is the only valid KOKKOS_ARCH. If no KOKKOS_ARCH is specified then"
  echo "*** 'BDW' will be used by default"
  echo "***"
  return
fi

export ATDM_CONFIG_SBATCH_DEFAULT_TIMEOUT=4:00:00

export ATDM_CONFIG_SLURM_DEFAULT_ACCOUNT=fy150090

source $ATDM_SCRIPT_DIR/common/toss3/environment_new.sh

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
