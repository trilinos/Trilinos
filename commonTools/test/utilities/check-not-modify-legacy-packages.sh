#!/bin/bash    

set -o nounset
set -o errexit
set -o pipefail

target_branch=${1}
source_branch=${2}

deprecated_packages=$(cat $(git rev-parse --show-toplevel)/CMakeLists.txt | grep "set(DEPRECATED_PACKAGES" | sed 's/set(DEPRECATED_PACKAGES //g' | sed 's/)//g' | tr '[:upper:]' '[:lower:]')
paths=
for package in ${deprecated_packages}
do
    paths="packages/${package} ${paths}"
done

estat=0
for commit in $(git log --format=%H ${source_branch} --not ${target_branch})
do
    echo "Processing commit ${commit}"
    git diff --quiet ${commit}^ ${commit} -- ${paths} || \
    {
        diff_output=$(git diff --name-only ${commit}^ ${commit} -- ${paths})
        for path in ${paths}
        do
            if [[ "${diff_output}" == *"${path}"* ]]
            then
                echo -e "  This commit modified a legacy package (${path}), which is deprecated and will be deleted soon."
            fi
        done
        estat=1
    }
done

exit ${estat}
