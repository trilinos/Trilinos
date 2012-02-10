#!/bin/bash

# Usage: setup-shared-repo.sh grpname

GROUP_NAME=$1
#echo $GROUP_NAME
OWNER_NAME=$2

chmod -R g+ws $PWD
chgrp -R $GROUP_NAME $PWD
if [ "$OWNER_NAME" != "" ] ; then
  chown -R $OWNER_NAME $PWD
fi
git repo-config core.sharedRepository true
