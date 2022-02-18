#!/usr/bin/env bash
SCRIPTFILE=$(realpath ${WORKSPACE:?}/Trilinos/packages/framework/pr_tools/PullRequestLinuxDriver.sh)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH:?}/common.bash
# set -x  # echo commands



# Load the right version of Git / Python based on a regex
# match to the Jenkins job name.
function bootstrap_modules() {
    print_banner "Bootstrap environment modules start"

    cuda_regex=".*(_cuda_).*"
    ride_regex=".*(ride).*"
    vortex_regex=".*(vortex).*"
    if [[ ${JOB_BASE_NAME:?} =~ ${cuda_regex} ]]; then
        if [[ ${NODE_NAME:?} =~ ${ride_regex} ]]; then
            message_std "PRDriver> " "Job is CUDA"
            execute_command_checked "module unload git"
            execute_command_checked "module unload python"
            execute_command_checked "module load git/2.10.1"
            execute_command_checked "module load python/3.7.3"
            get_python_packages pip3
            envvar_set_or_create PYTHON_EXE python3
        elif [[ ${NODE_NAME:?} =~ ${vortex_regex} ]]; then
            echo -e "Job is CUDA node is vortex"
            execute_command_checked "module load git/2.20.0"
            execute_command_checked "module load python/3.7.2"
            get_python_packages pip3
            envvar_set_or_create PYTHON_EXE python3
        else
            message_std "PRDriver> " "ERROR: Unable to find matching environment for CUDA job not on Ride."
            exit -1
        fi
    else
        execute_command_checked "module load apps/anaconda3.7"
        source /projects/sems/modulefiles/utils/sems-archive-modules-init.sh
        execute_command_checked "module unload sems-archive-git"
        execute_command_checked "module unload sems-archive-python"
        execute_command_checked "module load sems-archive-git/2.10.1"

        envvar_set_or_create     PYTHON_EXE $(which python3)
    fi

    module list

    message_std "PRDriver> " "Python EXE : ${PYTHON_EXE:?}"
    message_std "PRDriver> " "           : $(which ${PYTHON_EXE})"

    print_banner "Bootstrap environment modules complete"
}


print_banner "PullRequestLinuxDriver.sh"

# Set up Sandia PROXY environment vars
envvar_set_or_create https_proxy 'http://proxy.sandia.gov:80'
envvar_set_or_create http_proxy  'http://proxy.sandia.gov:80'
envvar_set_or_create no_proxy    'localhost,.sandia.gov,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'
#export https_proxy=http://proxy.sandia.gov:80
#export http_proxy=http://proxy.sandia.gov:80
#export no_proxy='localhost,.sandia.gov,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'


# bootstrap the python and git modules for this system
bootstrap_modules


# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH:?}/../..`
test -d ${REPO_ROOT:?}/.git || REPO_ROOT=`readlink -f ${WORKSPACE:?}/Trilinos`
message_std "PRDriver> " "REPO_ROOT : ${REPO_ROOT}"

# Get the md5 checksum of this script:
sig_script_old=$(get_md5sum ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriver.sh)

# Get the md5 checksum of the Merge script
sig_merge_old=$(get_md5sum ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriverMerge.py)


print_banner "Merge Source into Target"
message_std "PRDriver> " "TRILINOS_SOURCE_SHA: ${TRILINOS_SOURCE_SHA:?}"

# Prepare the command for the MERGE operation
merge_cmd_options=(
    ${TRILINOS_SOURCE_REPO:?}
    ${TRILINOS_SOURCE_BRANCH:?}
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
#err=$?
#if [ $err != 0 ]; then
#    print_banner "An error occurred during merge"
#    exit $err
#fi
print_banner "Merge completed"


print_banner "Check for PR Driver Script Modifications"

# Get the md5 checksum of this script:
#sig_script_new=$(get_md5sum ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriver.sh)
sig_script_new=$(get_md5sum ${SCRIPTFILE:?})
message_std "PRDriver> " ""
message_std "PRDriver> " "Script File: ${SCRIPTFILE:?}"
message_std "PRDriver> " "Old md5sum : ${sig_script_old:?}"
message_std "PRDriver> " "New md5sum : ${sig_script_new:?}"

# Get the md5 checksum of the Merge script
#sig_merge_new=$(get_md5sum ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriverMerge.py)
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
    ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriver.sh
    exit $?
fi

message_std "PRDriver> " ""
message_std "PRDriver> " "Driver and Merge scripts unchanged, proceeding to TEST phase"
message_std "PRDriver> " ""

# determine what MODE we are using
mode="standard"
if [[ "${JOB_BASE_NAME:?}" == "Trilinos_pullrequest_gcc_8.3.0_installation_testing" ]]; then
    mode="installation"
fi


envvar_set_or_create TRILINOS_BUILD_DIR ${WORKSPACE}/pull_request_test

#message_std "PRDriver> " "Create build directory if it does not exist."
#message_std "PRDriver> " "Build Dir: ${TRILINOS_BUILD_DIR:?}"
#mkdir -p ${TRILINOS_BUILD_DIR:?}



print_banner "Launch the Test Driver"


# Prepare the command for the TEST operation
test_cmd_options=(
    --source-repo-url=${TRILINOS_SOURCE_REPO:?}
    --source-branch-name=${TRILINOS_SOURCE_BRANCH:?}
    --target-repo-url=${TRILINOS_TARGET_REPO:?}
    --target-branch-name=${TRILINOS_TARGET_BRANCH:?}
    --pullrequest-build-name=${JOB_BASE_NAME:?}
    --genconfig-build-name=${GENCONFIG_BUILD_NAME:?}
    --pullrequest-env-config-file=${LOADENV_CONFIG_FILE:?}
    --pullrequest-gen-config-file=${GENCONFIG_CONFIG_FILE:?}
    --pullrequest-number=${PULLREQUESTNUM:?}
    --jenkins-job-number=${BUILD_NUMBER:?}
    --req-mem-per-core=3.0
    --max-cores-allowed=${TRILINOS_MAX_CORES:=29}
    --num-concurrent-tests=4
    --test-mode=${mode}
    --workspace-dir=${WORKSPACE:?}
    --filename-packageenables=${WORKSPACE:?}/packageEnables.cmake
    --filename-subprojects=${WORKSPACE:?}/package_subproject_list.cmake
    --source-dir=${WORKSPACE}/Trilinos
    --build-dir=${TRILINOS_BUILD_DIR:?}
    --ctest-driver=${WORKSPACE:?}/pr-ctest-framework/cmake/ctest-driver.cmake
    --ctest-drop-site=${TRILINOS_CTEST_DROP_SITE:?}
    #--dry-run
)



# Execute the TEST operation
test_cmd="${PYTHON_EXE:?} ${REPO_ROOT:?}/packages/framework/pr_tools/PullRequestLinuxDriverTest.py ${test_cmd_options[@]}"

# Call the script to launch the tests
print_banner "Execute Test Command"
message_std "PRDriver> " "cd $(pwd)"
message_std "PRDriver> " "${test_cmd:?} --pullrequest-cdash-track='${PULLREQUEST_CDASH_TRACK:?}'"
execute_command_checked "${test_cmd:?} --pullrequest-cdash-track='${PULLREQUEST_CDASH_TRACK:?}'"

#${test_cmd} --pullrequest-cdash-track="${PULLREQUEST_CDASH_TRACK:?}"
#exit $?
