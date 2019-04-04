ATDM_UTIL_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
unset ATDM_CONFIG_JOB_NAME
unset ATDM_CONFIG_BUILD_NAME
unset ATDM_CONFIG_SCRIPT_DIR
unset ATDM_CONFIG_TRILNOS_DIR
source ${ATDM_UTIL_SCRIPT_DIR}/unset_atdm_config_vars_system_name.sh
source ${ATDM_UTIL_SCRIPT_DIR}/unset_atdm_config_vars_build_options.sh
source ${ATDM_UTIL_SCRIPT_DIR}/unset_atdm_config_vars_environment.sh
