#!/bin/bash    

set -o nounset
set -o errexit
set -o pipefail

target_branch=${1}
source_branch=${2}

for commit in $(git log --format=%H ${source_branch} --not ${target_branch})
do
    echo "Processing commit ${commit}"
    git diff -U0 --ignore-all-space ${commit}~1 ${commit} | grep "^\+oid sha256" \
    && { echo "Git LFS pointers found; not allowed to add LFS files to Trilinos!" ; exit 1 ; }
done

echo "No LFS pointers found!"
