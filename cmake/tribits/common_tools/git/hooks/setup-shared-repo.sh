#!/bin/bash

# Usage: setup-shared-repo.sh grpname

GROUP_NAME=$1
#echo $GROUP_NAME
OWNER_NAME=$2
#echo $OWNER_NAME

echo "Making this a shared repo in git."
git repo-config core.sharedRepository true
echo "Setting the group sticky bit on directories so new files get the parent directory group."
find . -type d | xargs chmod g+s
echo "Removing the group sticky bit from files."
find . -type f | xargs chmod g-s
echo "Making repo group readable"
chmod -R g+r .
echo "Removing read/write/exec from others."
chmod -R o-rwxs .
echo "Changing group to $GROUP_NAME"
chgrp -R $GROUP_NAME .
if [ "$OWNER_NAME" != "" ] ; then
  echo "Changing owner on the git repo to to $OWNER_NAME"
  chown -R $OWNER_NAME .
  echo "Changing default owner on config to $OWNER_NAME"
  chown $OWNER_NAME config
fi
