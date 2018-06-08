################################################################################
#
# Set up env on serrano for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

# Chama and serrano jobs all use the same environmnet changes to the
# sourced script below will impact jobs on both of those
# machines. please be mindful of this when making changes

export ATDM_CONFIG_KOKKOS_ARCH=BDW
source $ATDM_SCRIPT_DIR/toss3/environment.sh
