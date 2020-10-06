#!/usr/bin/env bash
SCRIPTFILE=$(realpath $BASH_SOURCE)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH:?}/common.bash
# set -x  # echo commands


# TODO: Remove this
# Get the md5sum of a filename.
# param1: filename
# returns: md5sum of the file.
#function get_md5sum() {
#    local filename=${1:?}
#    local sig=$(md5sum ${filename:?} | cut -d' ' -f1)
#    echo "${sig:?}"
#}
#
# Get pip
# - @param1 python_exe - the python executable to install PIP for
#function get_pip() {
#    local python_exe=${1:?}
#
#    echo -e "--- Python: ${python_exe:?}"
#
#    # fetch get-pip.py
#    local curl_cmd="curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py"
#    echo -e "--- ${curl_cmd}"
#    eval ${curl_cmd}
#
#    get_pip_args=(
#        --user
#        --proxy="http://wwwproxy.sandia.gov:80"
#        --no-setuptools
#        --no-wheel
#    )
#    echo -e ""
#    echo -e "--- ${python_exe:?} ./get-pip.py ${get_pip_args[@]}"
#    ${python_exe:?} ./get-pip.py ${get_pip_args[@]}
#}



# Load the right version of Git / Python based on a regex
# match to the Jenkins job name.
function bootstrap_modules() {
    echo -e "PRDriver> ---------------------------------------"
    echo -e "PRDriver> Bootstrap environment modules Start"
    echo -e "PRDriver> ---------------------------------------"

    cuda_regex=".*(_cuda_).*"
    ride_regex=".*(ride).*"
    if [[ ${JOB_BASE_NAME:?} =~ ${cuda_regex} ]]; then
        if [[ ${NODE_NAME:?} =~ ${ride_regex} ]]; then
            echo -e "PRDriver> Job is CUDA"
            module unload git
            module unload python
            module load git/2.10.1
            #module load python/2.7.12
            module load python/3.7.3
            #get_pip python3
            get_python_packages pip3
            export PYTHON_EXE=python3
        else
            echo -e "PRDriver> ERROR: Unable to find matching environment for CUDA job not on Ride."
            exit -1
        fi
    else
        source /projects/sems/modulefiles/utils/sems-modules-init.sh
        module unload sems-git
        module unload sems-python
        module load sems-git/2.10.1
        module load sems-python/2.7.9
        # module load sems-python/3.5.2      # Currently not on cloud nodes
        #pip3 install --user configparser
        get_pip python2
        get_python_packages ${HOME}/.local/bin/pip2
        export PYTHON_EXE=python2
    fi

    module list

    echo -e "PRDriver> ---------------------------------------"
    echo -e "PRDriver> Bootstrap environment modules Complete"
    echo -e "PRDriver> ---------------------------------------"
}




echo -e "PRDRiver> ================================================="
echo -e "PRDriver> ="
echo -e "PRDriver> = PullRequestLinuxDriver.sh"
echo -e "PRDriver> ="
echo -e "PRDriver> ================================================="

# Set up Sandia PROXY environment vars
export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
export no_proxy='localhost,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'


# bootstrap the python and git modules for this system
bootstrap_modules


# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH:?}/../..`
echo -e "PRDriver> REPO_ROOT : ${REPO_ROOT}"

# Get the md5 checksum of this script:
sig_script_old=$(get_md5sum ${SCRIPTFILE:?})

# Get the md5 checksum of the Merge script
sig_merge_old=$(get_md5sum ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py)

# Prepare the command for the MERGE operation
merge_cmd_options=(
    ${TRILINOS_SOURCE_REPO:?}
    ${TRILINOS_SOURCE_BRANCH:?}
    ${TRILINOS_TARGET_REPO:?}
    ${TRILINOS_TARGET_BRANCH:?}
    ${TRILINOS_SOURCE_SHA:?}
    ${WORKSPACE:?}
    )
merge_cmd="python ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py ${merge_cmd_options[@]}"


# Call the script to handle merging the incoming branch into
# the current trilinos/develop branch for testing.
echo -e "PRDriver> "
echo -e "PRDriver> Execute Merge Command: ${merge_cmd:?}"
echo -e "PRDriver> "
${merge_cmd:?}
err=$?
if [ $err != 0 ]; then
    echo -e "PRDriver> An error occurred during merge"
    exit $err
else
    echo -e "PRDriver> Merge completed successfully."
fi
echo -e "PRDriver> "


# Get the md5 checksum of this script:
sig_script_new=$(get_md5sum ${SCRIPTFILE:?})
echo -e "PRDriver> Old md5 checksum ${sig_script_old:?} for ${SCRIPTFILE:?}"
echo -e "PRDriver> New md5 checksum ${sig_script_new:?} for ${SCRIPTFILE:?}"
echo -e "PRDriver> "

# Get the md5 checksum of the Merge script
sig_merge_new=$(get_md5sum ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py)
echo -e "PRDriver> Old md5 checksum ${sig_merge_old:?} for ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py"
echo -e "PRDriver> New md5 checksum ${sig_merge_new:?} for ${SCRIPTPATH}/PullRequestLinuxDriverMerge.py"

if [ "${sig_script_old:?}" != "${sig_script_new:?}" ] || [ "${sig_merge_old:?}" != "${sig_merge_new:?}"  ]
then
    echo -e "PRDriver> "
    echo -e "PRDriver> Driver or Merge script change detected. Re-launching PR Driver"
    echo -e "PRDriver> "
    ${SCRIPTFILE:?}
    exit $?
fi

echo -e "PRDriver> "
echo -e "PRDriver> Driver and Merge scripts unchanged, proceeding to TEST phase"
echo -e "PRDriver> "

# determine what MODE we are using
mode="standard"
if [[ "${JOB_BASE_NAME}" == "Trilinos_pullrequest_gcc_8.3.0_installation_testing" ]]; then
    mode="installation"
fi



# Prepare the command for the TEST operation
test_cmd_options=(
    --sourceRepo=${TRILINOS_SOURCE_REPO:?}
    --sourceBranch=${TRILINOS_SOURCE_BRANCH:?}
    --targetRepo=${TRILINOS_TARGET_REPO:?}
    --targetBranch=${TRILINOS_TARGET_BRANCH:?}
    --job_base_name=${JOB_BASE_NAME:?}
    --workspaceDir=${WORKSPACE:?}
    --github_pr_number=${PULLREQUESTNUM:?}
    --job_number=${BUILD_NUMBER:?}
    --req-mem-per-core=3.0
    --max-cores-allowed=29
    --num-concurrent-tests=4
    --mode=${mode}
    #--mode=installation
    --config="Trilinos/cmake/std/pr_config/pullrequest.ini"
    #--dry-run
)

# Execute the TEST operation
test_cmd="${PYTHON_EXE} ${SCRIPTPATH}/PullRequestLinuxDriverTest.py ${test_cmd_options[@]}"


# Call the script to launch the tests
echo -e "PRDriver> "
echo -e "PRDriver> Execute Test Command:"
echo -e "PRDriver> ${test_cmd:?}"
echo -e "PRDriver> "
${test_cmd}
exit $?



