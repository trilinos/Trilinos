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

cd $_SCRIPT_DIR
echo "Copy from: $PWD"

function copy_over_readonly {
  from_dir=$1
  file=$2
  to_dir=$3
  if [ -e $to_dir/$file ] ; then
    #echo "$to_dir/$file"
    chmod u+w $to_dir/$file
  fi
  cp -u -v $from_dir/$file $to_dir/$file
  file_lc=$(echo $file | tr '[:upper:]' '[:lower:]')
  if [ -e $to_dir/$file_lc ] ; then
    rm -f $to_dir/$file_lc
  fi
  cp -u -v $from_dir/$file $to_dir/$file_lc

}


function create_symlink {
  file_from=$1
  file_to=$2
  if [ -e $file_to ] ; then
    echo "rm -f $file_to"
    rm -f $file_to
  fi
  echo "ln -s $file_from $file_to"
  ln -s $file_from $file_to
}


function rm_if_exists {
  file_to_rm=$1
  if [ -e $file_to_rm ] ; then
    echo "Removing old file $file_to_rm"
    rm -f $file_to_rm
  fi

}



_DEST_BASE_DIR=$1

echo "Copy to: $_DEST_BASE_DIR" 
copy_over_readonly  guides/users_guide  TribitsUsersGuide.html  $_DEST_BASE_DIR
copy_over_readonly  guides/maintainers_guide  TribitsMaintainersGuide.html  $_DEST_BASE_DIR
copy_over_readonly  build_ref  TribitsBuildReference.html  $_DEST_BASE_DIR
copy_over_readonly  build_ref  TribitsBuildReference.pdf  $_DEST_BASE_DIR

echo "Create symlinks for old TribitsDevelopersGuide.html"
cd $_DEST_BASE_DIR
create_symlink  TribitsUsersGuide.html  TribitsDevelopersGuide.html
create_symlink  TribitsUsersGuide.html  tribitsdevelopersguide.html

echo "Removing old copies for PDf files for now"
rm_if_exists TribitsDevelopersGuide.pdf
rm_if_exists tribitsdevelopersguide.pdf
