#!/bin/bash -e
#
# Run this script to set up a shared set of directories/files as:
#
#  <this-dir>/setup-shared-dir.sh <base-dir> <owning-group>
#
# This script set up permissions to allow anyone in <owning-group> to
# read-write the files and directories under <base-dir> to allow them to
# maintain the files.  The permissions on the files/directories are to allow
# anyone on system to read (or execute) the files, but not modify them.
#
# NOTE: If you must limit who can access the files under <base-dir>, then you
# can't used this script or if you do, you will need to lock down the
# permissions after you run it.
#
# NOTE: If you don't have permissions to write all of these files and
# directories you might need to run this script as:
#
#  sudo <this-dir>/setup-shared-dir.sh <base-dir> <owning-group>
#

BASE_DIR=$1
OWNING_GROUP=$2
echo "Standardizing the ownership and permissions of $BASE_DIR with owning group $OWNING_GROUP"

# Set the owning group for everything
chgrp -R $OWNING_GROUP $BASE_DIR
# Give the owning group read/write for all files
chmod -R g+rw $BASE_DIR
# Allow everyone to read everything and execute if it is executable!
chmod -R a+rX $BASE_DIR
