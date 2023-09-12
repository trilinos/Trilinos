#!/bin/bash    

set -o nounset
set -o errexit
set -o pipefail

claimed_start_of_branch=${1}
end_of_branch=${2}
start_of_branch=$(git merge-base ${claimed_start_of_branch} ${end_of_branch})

for commit in $(git log --format=%H ${start_of_branch}..${end_of_branch})
do
    echo "Processing commit ${commit}"
    git diff -U0 --ignore-all-space ${commit}~1 ${commit} | grep "^\+oid sha256" \
    && { echo "Git LFS pointers found; not allowed to add LFS files to Trilinos!" ; exit 1 ; }
done

echo "No LFS pointers found!"
