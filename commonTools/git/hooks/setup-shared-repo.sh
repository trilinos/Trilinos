#!/bin/bash

# Usage: setup-shared-repo.sh grpname

GROUP_NAME=$1
#echo $GROUP_NAME

chmod -R g+ws *
chgrp -R $GROUP_NAME *
git repo-config core.sharedRepository true
