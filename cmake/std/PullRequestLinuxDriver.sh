#!/usr/bin/env bash
# set -x  # echo commands

function get_scriptname() {
    # Get the full path to the current script
    local script_name=`basename $0`
    local script_path=$(dirname $(readlink -f $0))
    local script_file="${script_path}/${script_name:?}"
    echo "${script_file}"
}

function get_scriptpath() {
    # Get the full path to the current script
    local script_name=`basename $0`
    local script_path=$(dirname $(readlink -f $0))
    echo "${script_path}"
}

# Get the md5sum of a filename.
# param1: filename
# returns: md5sum of the file.
function get_md5sum() {
    local filename=${1:?}
    local sig=$(md5sum ${filename:?} | cut -d' ' -f1)
    echo "${sig:?}"
}


# Set up Sandia PROXY environment vars
export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
export no_proxy='localhost,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'


# Load the right version of Git / Python based on a regex 
# match to the Jenkins job name.
cuda_regex=".*(_cuda_).*"
ride_regex=".*(ride).*"
if [[ ${JOB_BASE_NAME:?} =~ ${cuda_regex} ]]; then
    if [[ ${NODE_NAME:?} =~ ${ride_regex} ]]; then
        echo -e "Job is CUDA"
        module load git/2.10.1
        module load python/2.7.12
    else
        echo -e "ERROR: Unable to find matching environment for CUDA job not on Ride."
        exit -1
    fi
else
    source /projects/sems/modulefiles/utils/sems-modules-init.sh
    module load sems-git/2.10.1
    module load sems-python/2.7.9
fi


# Identify the path to this script
SCRIPTPATH=$(get_scriptpath)
script_file=$(get_scriptname)


# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH}/../..`
echo -e "REPO_ROOT : ${REPO_ROOT}"

# Get the md5 checksum of this script:
sig_script_old=$(get_md5sum ${script_file:?})
echo -e ""
echo -e ">>> Old md5 checksum for '${script_file:?}' = ${sig_script_old}"
echo -e ""

# Get the md5 checksum of the Merge script
sig_merge_old=$(get_md5sum ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py)
echo -e ""
echo -e ">>> Old md5 checksum for '${SCRIPTPATH}/PullRequestLinuxDriverMerge.py' = ${sig_merge_old}"
echo -e ""

# Prepare the command for the MERGE operation
merge_cmd_options=(
    ${TRILINOS_SOURCE_REPO:?}
    ${TRILINOS_SOURCE_BRANCH:?}
    ${TRILINOS_TARGET_REPO:?}
    ${TRILINOS_TARGET_BRANCH:?}
    ${TRILINOS_SOURCE_SHA:?}
    ${WORKSPACE:?}
    )
merge_cmd=${SCRIPTPATH}/PullRequestLinuxDriverMerge.py ${merge_cmd_options[@]}


# Call the script to handle merging the incoming branch into
# the current trilinos/develop branch for testing.
echo -e ""
echo -e "Execute Merge Command:"
echo -e "${merge_cmd:?} 
echo -e ""
${merge_cmd:?}


# Get the md5 checksum of this script:
sig_script_new=$(get_md5sum ${script_file:?})
echo ""
echo -e ">>> New md5 checksum for '${script_file:?}' = ${sig_script_new}"
echo ""

# Get the md5 checksum of the Merge script
sig_merge_new=$(get_md5sum ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py)
echo ""
echo -e ">>> New md5 checksum for '${SCRIPTPATH}/PullRequestLinuxDriverMerge.py' = ${sig_merge_new}"
echo ""

if [ "${sig_script_old:?}" != "${sig_script_new:?}" ] || [ "${sig_merge_old:?}" != "${sig_merge_new:?}"  ]
then
    echo -e ""
    echo "Driver or Merge script change detected. Re-launching PR Driver"
    echo -e ""
    ${merge_cmd:?}
    exit $?
fi

echo -e ""
echo -e "Driver and Merge scripts unchaged, proceeding to TEST phase"
echo -e ""

# Prepare the command for the TEST operation
test_cmd_options=(
    ${TRILINOS_SOURCE_REPO:?}
    ${TRILINOS_SOURCE_BRANCH:?}
    ${TRILINOS_TARGET_REPO:?}
    ${TRILINOS_TARGET_BRANCH:?}
    ${JOB_BASE_NAME:?}
    ${PULLREQUESTNUM:?}
    ${BUILD_NUMBER:?}
    ${WORKSPACE:?}
    )
test_cmd=${SCRIPTPATH}/PullRequestLinuxDriverTest.py ${test_cmd_options[@]}

# Call the script to launch the tests
echo -e ""
echo -e "Execute Test Command:"
echo -e "${test_cmd:?} 
echo -e ""
${test_cmd}


