#!/bin/bash

# Usage: setup-shared-repo.sh grpname

GROUP_NAME=$1
#echo $GROUP_NAME

chmod -R g+ws $PWD
chgrp -R $GROUP_NAME $PWD
git repo-config core.sharedRepository true
