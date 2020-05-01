# These are defaults that work for almost every system where that ATDM DevOps
# effort installs Trilinos as a biproduct of that ATDM Trilinos builds that
# submit to CDash.  This assumes that either the 'jenkins' or the
# 'atdm-devops-admin' entity accounts is driving the ctest -S scripts and the
# group 'wg-run-as-atdm-devsops' exists on the machine.

if [[ "${ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT}" == "" ]] ; then
  export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=/projects/atdm_devops/trilinos_installs
fi

if [[ "${ATDM_CONFIG_MAKE_INSTALL_GROUP_DEFAULT}" == "" ]] ; then
  export ATDM_CONFIG_MAKE_INSTALL_GROUP_DEFAULT=wg-run-as-atdm-devops
fi
