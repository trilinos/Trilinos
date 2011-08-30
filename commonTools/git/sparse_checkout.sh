#!/bin/bash

# Run this script in the main Trilinos repo dir to do a sparse
# checkout of directories.  For example:
#
#  $ cd Trilinos
#  $ ./commonTools/git/sparse_checkout.sh commonTools packages
#  $ ls
#
# should return just 'commonTools' and 'packages'
#
# Warning: Make sure you run this script before you change any files
# or the 'git read-tree ...' command will blow them away (a lesson I
# learned the hard way).  In general, commit your work *before* your
# run this script.

DIRS_FILES_LIST=$@

# Tell Git to allow sparse checkout
git config core.sparsecheckout true

# Set up the list of dirs/files for sparse checkout
SC_FILE=.git/info/sparse-checkout
echo > $SC_FILE
for dir_or_file in $DIRS_FILES_LIST ; do
  echo "File/directory to keep: $dir_or_file";
  echo $dir_or_file >> $SC_FILE
done

# Tell git to process the list of dirs/files
git read-tree -m -u HEAD
