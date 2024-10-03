#!/usr/bin/env bash
SCRIPTFILE=$(realpath ${WORKSPACE:?}/Trilinos/packages/framework/pr_tools/PullRequestLinuxDriver.sh)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH:?}/common.bash


# Configure ccache via environment variables
function configure_ccache() {
    print_banner "Configuring ccache"

    envvar_set_or_create CCACHE_NODISABLE true
    envvar_set_or_create CCACHE_DIR '/fgs/trilinos/ccache/cache'
    envvar_set_or_create CCACHE_BASEDIR "${WORKSPACE:?}"
    envvar_set_or_create CCACHE_NOHARDLINK true
    envvar_set_or_create CCACHE_UMASK 077
    envvar_set_or_create CCACHE_MAXSIZE 100G

    message_std "PRDriver> " "$(ccache --show-stats --verbose)"
}


# Load the right version of Git / Python based on a regex
# match to the Jenkins job name.
function bootstrap_modules() {
    print_banner "Bootstrap environment modules start"
    message_std "PRDriver> " "Job is $JOB_BASE_NAME"

    vortex_regex=".*(vortex).*"
    container_regex=".*(container).*"
    if [[ ${NODE_NAME:?} =~ ${vortex_regex} || ${on_ats2} == "1" ]]; then
        execute_command_checked "module load git/2.20.0"
        execute_command_checked "module load python/3.7.2"
        get_python_packages pip3
        # Always create user's tmp dir for nvcc. See https://github.com/trilinos/Trilinos/issues/10428#issuecomment-1109956415.
        mkdir -p /tmp/trilinos

        module list
    elif [[ ${NODE_NAME:?} =~ ${container_regex} ]]; then
	echo "Nothing done for bootstrap in a container"
	module list
    elif [[ ${on_weaver} == "1" ]]; then
        module unload git
        module unload python
        module load git/2.10.1
        module load python/3.7.3
        get_python_packages pip3

        module list
    elif [[ ${on_rhel8} == "1" ]]; then
        source /projects/sems/modulefiles/utils/sems-modules-init.sh
        module unload sems-git
        module unload sems-python
        module load sems-git/2.37.0
        module load sems-python/3.9.0

        module list
    else
        source /projects/sems/modulefiles/utils/sems-modules-init.sh
        execute_command_checked "module unload sems-git"
        execute_command_checked "module unload sems-python"
        execute_command_checked "module load sems-git/2.37.0"
        execute_command_checked "module load sems-python/3.9.0"
        execute_command_checked "module load sems-ccache"
        configure_ccache

        module list
    fi

    print_banner "Bootstrap environment modules complete"
}


print_banner "PullRequestLinuxDriver.sh"

# Argument defaults
on_weaver=0
on_ats2=0
on_kokkos_develop=0
on_rhel8=0
bootstrap=1

original_args=$@

# Do POSIXLY_CORRECT option handling.
ARGS=$(getopt -n PullRequestLinuxDriver.sh \
 --options '+x' \
 --longoptions on-rhel8,on_rhel8 \
 --longoptions on-weaver,on_weaver \
 --longoptions on-ats2,on_ats2 \
 --longoptions kokkos-develop \
 --longoptions extra-configure-args: \
 --longoptions no-bootstrap -- "${@}") || exit $?

eval set -- "${ARGS}"

while [ "$#" -gt 0 ]
do
    case "${1}" in
    (--on_weaver|--on-weaver)
        on_weaver=1
        shift
        ;;
    (--on_rhel8|--on-rhel8)
        on_rhel8=1
        shift
        ;;
    (--on_ats2|--on-ats2)
        on_ats2=1
        shift
        ;;
    (--kokkos-develop)
        on_kokkos_develop=1
        shift
        ;;
    (--no-bootstrap)
        bootstrap=0
        shift
        ;;
    (--extra-configure-args)
        extra_configure_args=$2
        shift 2
        ;;
    (-h|--help)
        # When help is requested echo it to stdout.
        echo -e "$USAGE"
        exit 0
        ;;
    (-x)
        set -x
        shift
        ;;
    (--) # This is an explicit directive to stop processing options.
        shift
        break
        ;;
    (-*) # Catch options which are defined but not implemented.
        echo >&2 "${toolName}: ${1}: Unimplemented option passed."
        exit 1
        ;;
    (*) # The first parameter terminates option processing.
        break
        ;;
    esac
done

# Set up Sandia PROXY environment vars
envvar_set_or_create https_proxy 'http://proxy.sandia.gov:80'
envvar_set_or_create http_proxy  'http://proxy.sandia.gov:80'
envvar_set_or_create no_proxy    'localhost,.sandia.gov,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'

# bootstrap the python and git modules for this system
if [[ ${bootstrap} == "1" ]]; then
    bootstrap_modules
fi

envvar_set_or_create PYTHON_EXE $(which python3)
message_std "PRDriver> " "Python EXE : ${PYTHON_EXE:?}"

# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH:?}/../..`
test -d ${REPO_ROOT:?}/.git || REPO_ROOT=`readlink -f ${WORKSPACE:?}/Trilinos`
message_std "PRDriver> " "REPO_ROOT : ${REPO_ROOT}"

# Get the md5 checksum of this script:
sig_script_old=$(get_md5sum ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriver.sh)

# Get the md5 checksum of the Merge script
sig_merge_old=$(get_md5sum ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriverMerge.py)

if [[ ${on_kokkos_develop} == "1" ]]; then
    message_std "PRDriver> --kokkos-develop is set - setting kokkos and kokkos-kernels packages to current develop and pointing at them"
    "${SCRIPTPATH}"/SetKokkosDevelop.sh
    extra_configure_args="\"-DKokkos_SOURCE_DIR_OVERRIDE:string=kokkos;-DKokkosKernels_SOURCE_DIR_OVERRIDE:string=kokkos-kernels\"${extra_configure_args:+;${extra_configure_args}}"
else
    print_banner "Merge Source into Target"
    message_std "PRDriver> " "TRILINOS_SOURCE_SHA: ${TRILINOS_SOURCE_SHA:?}"

    # Prepare the command for the MERGE operation
    merge_cmd_options=(
        ${TRILINOS_SOURCE_REPO:?}
        ${TRILINOS_TARGET_REPO:?}
        ${TRILINOS_TARGET_BRANCH:?}
        ${TRILINOS_SOURCE_SHA:?}
        ${WORKSPACE:?}
        )
    merge_cmd="${PYTHON_EXE:?} ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriverMerge.py ${merge_cmd_options[@]}"


    # Call the script to handle merging the incoming branch into
    # the current trilinos/develop branch for testing.
    message_std "PRDriver> " ""
    message_std "PRDriver> " "Execute Merge Command: ${merge_cmd:?}"
    message_std "PRDriver> " ""
    execute_command_checked "${merge_cmd:?}"
    print_banner "Merge completed"


    print_banner "Check for PR Driver Script Modifications"

    # Get the md5 checksum of this script:
    sig_script_new=$(get_md5sum ${SCRIPTFILE:?})
    message_std "PRDriver> " ""
    message_std "PRDriver> " "Script File: ${SCRIPTFILE:?}"
    message_std "PRDriver> " "Old md5sum : ${sig_script_old:?}"
    message_std "PRDriver> " "New md5sum : ${sig_script_new:?}"

    # Get the md5 checksum of the Merge script
    export MERGE_SCRIPT=${SCRIPTPATH:?}/PullRequestLinuxDriverMerge.py
    sig_merge_new=$(get_md5sum ${MERGE_SCRIPT:?})
    message_std "PRDriver> " ""
    message_std "PRDriver> " "Script File: ${MERGE_SCRIPT:?}"
    message_std "PRDriver> " "Old md5sum : ${sig_merge_old:?}"
    message_std "PRDriver> " "New md5sum : ${sig_merge_new:?}"

    if [ "${sig_script_old:?}" != "${sig_script_new:?}" ] || [ "${sig_merge_old:?}" != "${sig_merge_new:?}"  ]
    then
        message_std "PRDriver> " ""
        message_std "PRDriver> " "Driver or Merge script change detected. Re-launching PR Driver"
        message_std "PRDriver> " ""
        ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriver.sh $original_args
        exit $?
    fi

    message_std "PRDriver> " ""
    message_std "PRDriver> " "Driver and Merge scripts unchanged, proceeding to TEST phase"
    message_std "PRDriver> " ""
fi

# determine what MODE we are using
mode="standard"
if [[ "${JOB_BASE_NAME:?}" == "Trilinos_pullrequest_gcc_8.3.0_installation_testing" ]]; then
    mode="installation"
fi

envvar_set_or_create TRILINOS_BUILD_DIR ${WORKSPACE}/pull_request_test

print_banner "Launch the Test Driver"

# Prepare the command for the TEST operation
test_cmd_options=(
    --target-branch-name=${TRILINOS_TARGET_BRANCH:?}
    --genconfig-build-name=${GENCONFIG_BUILD_NAME:?}
    --pullrequest-env-config-file=${LOADENV_CONFIG_FILE:?}
    --pullrequest-gen-config-file=${GENCONFIG_CONFIG_FILE:?}
    --pullrequest-number=${PULLREQUESTNUM:?}
    --jenkins-job-number=${BUILD_NUMBER:?}
    --req-mem-per-core=4.0
    --max-cores-allowed=${TRILINOS_MAX_CORES:=29}
    --num-concurrent-tests=16
    --test-mode=${mode}
    --workspace-dir=${WORKSPACE:?}
    --filename-packageenables=${WORKSPACE:?}/packageEnables.cmake
    --filename-subprojects=${WORKSPACE:?}/package_subproject_list.cmake
    --source-dir=${WORKSPACE}/Trilinos
    --build-dir=${TRILINOS_BUILD_DIR:?}
    --ctest-driver=${WORKSPACE:?}/Trilinos/cmake/SimpleTesting/cmake/ctest-driver.cmake
    --ctest-drop-site=${TRILINOS_CTEST_DROP_SITE:?}
    --dashboard-build-name=${DASHBOARD_BUILD_NAME}
)

if [[ ${extra_configure_args} ]]
then
    test_cmd_options+=( "--extra-configure-args=\"${extra_configure_args}\"")
fi

if [[ ${GENCONFIG_BUILD_NAME} == *"gnu"* ]]
then
    test_cmd_options+=( "--use-explicit-cachefile ")
fi

test_cmd="${PYTHON_EXE:?} ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriverTest.py ${test_cmd_options[@]}"

# Call the script to launch the tests
print_banner "Execute Test Command"
message_std "PRDriver> " "cd $(pwd)"
message_std "PRDriver> " "${test_cmd:?} --pullrequest-cdash-track='${PULLREQUEST_CDASH_TRACK:?}'"
execute_command_checked "${test_cmd:?} --pullrequest-cdash-track='${PULLREQUEST_CDASH_TRACK:?}'"
