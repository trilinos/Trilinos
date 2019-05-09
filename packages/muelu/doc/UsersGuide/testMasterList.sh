#!/bin/bash

diff_command="diff -w"

# copy all *.tex to *.tex.gold
shopt -s nullglob
for i in *.tex ; do
    cp $i $i.repo
done

# run update script, write MueLu_MasterList.cpp to this folder
./update_params.sh MueLu_MasterList.cpp

# check whether xsltproc crashed
if [ $? -ne 0 ]; then
    echo "xsltproc is available, but crashed when the update_params.sh script was run. Consider updating xsltproc."
    exit 0
fi

# for simplicity, compare all tex files
test_files="*.tex"
test_files+=" MueLu_MasterList.cpp"

return_code=0
for i in $test_files; do
    $diff_command $i $i.repo
    if [ $? -eq 0 ]; then
        echo "$i matches"
    else
        echo "$i does not match"
        return_code=1
    fi
done

if [ $return_code -eq 1 ]; then
    echo "Do not edit MueLu_MasterList.cpp directly, but edit masterList.xml and run update_params.sh"
fi

exit $return_code
