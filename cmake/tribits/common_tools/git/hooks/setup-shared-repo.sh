#!/bin/bash

# Usage: setup-shared-repo.sh grpname

GROUP_NAME=$1
#echo $GROUP_NAME

echo "Setting the group sticky bit on directories so new files get the parent directory group."
find . -type d | xargs chmod g+s
echo "Removing the group sticky bit from files."
find . -type f | xargs chmod g-s
if [ "$GROUP_NAME" == "WORLD" ] ; then
  echo "Tell git to make the shared git repo world read/write for new git files."
  git repo-config core.sharedRepository world
  echo "Setting all directories to be read/write/execute for world."
  find . -type d | xargs chmod a+rwx
  echo "Setting all files to be read/write for world."
  find . -type f | xargs chmod a+rw
else
  echo "Tell git to make the shared git repo group read/write for new git files."
  git repo-config core.sharedRepository true
  echo "Making repo group read/write"
  chmod -R g+rw .
  echo "Removing read/write/exec from others."
  chmod -R o-rwxs .
  echo "Changing group to $GROUP_NAME"
  chgrp -R $GROUP_NAME .
fi
