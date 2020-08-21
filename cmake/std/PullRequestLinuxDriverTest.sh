#!/usr/bin/env bash
# set -x  # echo commands

# Identify the path to this script
SCRIPTPATH="$(cd "$(dirname "$0")" ; pwd -P)"
echo -e "SCRIPTPATH: ${SCRIPTPATH}"

# Identify the path to the trilinos repository root
REPO_ROOT=`readlink -f ${SCRIPTPATH}/../..`
echo -e "REPO_ROOT : ${REPO_ROOT}"

# Both scripts will need access through the sandia
# proxies so set them here.
export https_proxy=http://wwwproxy.sandia.gov:80
export http_proxy=http://wwwproxy.sandia.gov:80
export no_proxy='localhost,localnets,127.0.0.1,169.254.0.0/16,forge.sandia.gov'


options=(
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
    --mode=installation
    --config="Trilinos/cmake/std/configs/trilinos_pr.ini"
    #--dry-run
)    

echo "=========================================================="
echo "    LAUNCH PullRequestLinuxDriverTest.py "
echo "which python3: $(which python3)"
echo "which pip3   : $(which pip3)"
echo "cmd:"
echo "python3 ${SCRIPTPATH}/PullRequestLinuxDriverTest.py ${options[@]}"

python3 ${SCRIPTPATH}/PullRequestLinuxDriverTest.py ${options[@]}
err=$?
if [ $err != 0 ]; then
    exit $err
fi

