#!/usr/bin/env bash
SCRIPTFILE=$(realpath ${WORKSPACE:?}/Trilinos/cmake/std/PullRequestLinuxDriver.sh)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH:?}/common.bash
# set -x  # echo commands



# Load the right version of Git / Python based on a regex
# match to the Jenkins job name.
function bootstrap_modules() {
    print_banner "Bootstrap environment modules start"

    cuda_regex=".*(_cuda_).*"
    weaver_regex=".*(weaver).*"
    vortex_regex=".*(vortex).*"
    if [[ ${JOB_BASE_NAME:?} =~ ${cuda_regex} ]]; then
        if [[ ${NODE_NAME:?} =~ ${weaver_regex} ]]; then
            message_std "PRDriver> " "Job is CUDA"
            module unload git
            module unload python
            module load git/2.10.1
            module load python/3.7.3
            get_python_packages pip3
            export PYTHON_EXE=python3
        elif [[ ${NODE_NAME:?} =~ ${vortex_regex} ]]; then
            echo -e "Job is CUDA node is vortex"
            module load git/2.20.0
            module load python/3.7.2
            get_python_packages pip3
            export PYTHON_EXE=python3
        else
            message_std "PRDriver> " "ERROR: Unable to find matching environment for CUDA job not on Ride."
            exit -1
        fi
    else
        source /projects/sems/modulefiles/utils/sems-archive-modules-init.sh
        module unload sems-archive-git
        module unload sems-archive-python
        module load sems-archive-git/2.10.1

#        module load sems-python/3.5.2      # Currently not on cloud nodes
#        #pip3 install --user configparser
#        get_python_packages pip3
#        export PYTHON_EXE=python3

#         envvar_set_or_create     PYTHONHOME /projects/sierra/linux_rh7/install/Python/3.6.3
         envvar_set_or_create     PYTHONHOME /projects/sierra/linux_rh7/install/Python/3.6.10
#        #envvar_set_or_create     PYTHONPATH ${HOME}/.local/lib/python3.6/site-packages
#        #envvar_append_or_create  PYTHONPATH ${PYTHONHOME:?}/lib/python3.6/site-packages
#        unset PYTHONHOME
         unset PYTHONPATH
         envvar_prepend_or_create PATH       ${PYTHONHOME:?}/bin
         envvar_set_or_create     PYTHON_EXE ${PYTHONHOME:?}/bin/python3
         #export PYTHONHOME=/projects/sierra/linux_rh7/install/Python/3.6.3
         #export PYTHONPATH=${HOME}/.local/lib/python3.6/site-packages:${PYTHONHOME:?}/lib/python3.6/site-packages
    fi

    module list

    message_std "PRDriver> " "Python EXE : ${PYTHON_EXE:?}"
    message_std "PRDriver> " "           : $(which ${PYTHON_EXE})"

    print_banner "Bootstrap environment modules complete"
}


print_banner "PullRequestLinuxDriver.sh"

# Set up Sandia PROXY environment vars
if [[ "${TRILINOS_PR_DO_NOT_SET_PROXY}}" == "" ]] ; then
  export https_proxy=http://proxy.sandia.gov:80
  export http_proxy=http://proxy.sandia.gov:80
  export no_proxy='localhost,.sandia.gov,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'
fi


# bootstrap the python and git modules for this system
bootstrap_modules


# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH:?}/../..`
test -d ${REPO_ROOT:?}/.git || REPO_ROOT=`readlink -f ${WORKSPACE:?}/Trilinos`
message_std "PRDriver> " "REPO_ROOT : ${REPO_ROOT}"

# Get the md5 checksum of this script:
sig_script_old=$(get_md5sum ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriver.sh)

# Get the md5 checksum of the Merge script
sig_merge_old=$(get_md5sum ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriverMerge.py)


print_banner "Merge Source into Target"

# Prepare the command for the MERGE operation
merge_cmd_options=(
    ${TRILINOS_SOURCE_REPO:?}
    ${TRILINOS_SOURCE_BRANCH:?}
    ${TRILINOS_TARGET_REPO:?}
    ${TRILINOS_TARGET_BRANCH:?}
    ${TRILINOS_SOURCE_SHA:?}
    ${WORKSPACE:?}
    )
merge_cmd="${PYTHON_EXE:?} ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriverMerge.py ${merge_cmd_options[@]}"


# Call the script to handle merging the incoming branch into
# the current trilinos/develop branch for testing.
message_std "PRDriver> " ""
message_std "PRDriver> " "Execute Merge Command: ${merge_cmd:?}"
message_std "PRDriver> " ""
${merge_cmd:?}
err=$?
if [ $err != 0 ]; then
    print_banner "An error occurred during merge"
    exit $err
fi
print_banner "Merge completed"


# Get the md5 checksum of this script:
sig_script_new=$(get_md5sum ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriver.sh)
message_std "PRDriver> " "Old md5 checksum ${sig_script_old:?} for ${SCRIPTFILE:?}"
message_std "PRDriver> " "New md5 checksum ${sig_script_new:?} for ${SCRIPTFILE:?}"
message_std "PRDriver> " ""

# Get the md5 checksum of the Merge script
sig_merge_new=$(get_md5sum ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriverMerge.py)
message_std "PRDriver> " "Old md5 checksum ${sig_merge_old:?} for ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py"
message_std "PRDriver> " "New md5 checksum ${sig_merge_new:?} for ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py"

if [ "${sig_script_old:?}" != "${sig_script_new:?}" ] || [ "${sig_merge_old:?}" != "${sig_merge_new:?}"  ]
then
    message_std "PRDriver> " ""
    message_std "PRDriver> " "Driver or Merge script change detected. Re-launching PR Driver"
    message_std "PRDriver> " ""
    ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriver.sh
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


print_banner "Launch the Test Driver"

# Prepare the command for the TEST operation
test_cmd_options=(
    --source-repo-url=${TRILINOS_SOURCE_REPO:?}
    --source-branch-name=${TRILINOS_SOURCE_BRANCH:?}
    --target-repo-url=${TRILINOS_TARGET_REPO:?}
    --target-branch-name=${TRILINOS_TARGET_BRANCH:?}
    --pullrequest-build-name=${JOB_BASE_NAME:?}
    --pullrequest-config-file="Trilinos/cmake/std/pr_config/pullrequest.ini"
    --pullrequest-number=${PULLREQUESTNUM:?}
    --jenkins-job-number=${BUILD_NUMBER:?}
    --req-mem-per-core=3.0
    --max-cores-allowed=29
    --num-concurrent-tests=4
    --test-mode=${mode}
    --workspace-dir=${WORKSPACE:?}
    #--dry-run
)

# Execute the TEST operation
test_cmd="${PYTHON_EXE:?} ${REPO_ROOT:?}/cmake/std/PullRequestLinuxDriverTest.py ${test_cmd_options[@]}"


# Call the script to launch the tests
print_banner "Execute Test Command"
message_std "PRDriver> " "cd $(pwd)"
message_std "PRDriver> " "${test_cmd:?} --pullrequest-cdash-track='${PULLREQUEST_CDASH_TRACK:?}'"
${test_cmd} --pullrequest-cdash-track="${PULLREQUEST_CDASH_TRACK:?}"
exit $?
