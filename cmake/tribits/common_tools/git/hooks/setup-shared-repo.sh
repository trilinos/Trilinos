#!/bin/bash

# Usage: setup-shared-repo.sh grpname

GROUP_NAME=$1
#echo $GROUP_NAME
OWNER_NAME=$2
#echo $OWNER_NAME

echo "Making this a shared repo in git."
git repo-config core.sharedRepository true
echo "Making repo group read/writable and sticky" 
chmod -R g+rws .
echo "Changing group to $GROUP_NAME"
chgrp -R $GROUP_NAME .
if [ "$OWNER_NAME" != "" ] ; then
  echo "Changing default owner on . to $OWNER_NAME"
  chown -R $OWNER_NAME .
  echo "Changing default owner on config to $OWNER_NAME"
  chown $OWNER_NAME config
fi
