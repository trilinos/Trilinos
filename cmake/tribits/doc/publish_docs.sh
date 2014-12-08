#!/bin/bash

# This simple script copies selected TriBITS documents (built using the script
# 'build_docs.sh') so they can be viewed through a web browser.  This is run
# as:
#
#    <some-base-dir>/publish_docs.sh <destination-dir>
#
# (where <destination-dir> must be an absolute path).  The script only needs
# to be given the (already existing) destination base directory as an
# argument.  The script already knows the location of the TriBITS source
# directories if this is run out of the source tree.

#_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
_ABS_FILE_PATH=`readlink -f $0`
#echo "_ABS_FILE_PATH = '$_ABS_FILE_PATH'"
_SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"
_TRIBITS_BASE_DIR=$_SCRIPT_DIR/..
echo "Copy from: $_TRIBITS_BASE_DIR"

function copy_over_readonly {
  from_dir=$1
  file=$2
  to_dir=$3
  if [ -e $to_dir/$file ] ; then
    #echo "$to_dir/$file"
    chmod u+w $to_dir/$file
  fi
  cp -u -v $from_dir/$file $to_dir/$file
}

_DEST_BASE_DIR=$1
echo "Copy to: $_DEST_BASE_DIR" 
copy_over_readonly  developers_guide  TribitsDevelopersGuide.html  $_DEST_BASE_DIR
copy_over_readonly  developers_guide  TribitsDevelopersGuide.pdf  $_DEST_BASE_DIR
copy_over_readonly  build_quick_ref  TribitsBuildQuickRef.html  $_DEST_BASE_DIR
copy_over_readonly  build_quick_ref  TribitsBuildQuickRef.pdf  $_DEST_BASE_DIR
