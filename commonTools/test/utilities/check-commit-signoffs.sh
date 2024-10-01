#!/bin/bash    

set -o nounset
set -o errexit
set -o pipefail

target_branch=${1}
source_branch=${2}

estat=0
for commit in $(git log --format=%H ${source_branch} --not ${target_branch})
do
    echo "Processing commit ${commit}"
    git show -s --format=%B ${commit} | grep --extended-regexp --quiet "Signed-off-by:\s+\S+.*<\S+@\S+\.\S+>" \
    || { echo -e "Commit ${commit} does not contain the required DCO (https://developercertificate.org) sign-off!\nThe \"DCO\" check for this PR should have failed, and manual override is not permitted.\n" ; estat=1 ; }
done

exit ${estat}
