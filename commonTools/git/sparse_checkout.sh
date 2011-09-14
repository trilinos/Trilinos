#!/bin/bash

# Run this script in the main Trilinos repo dir to do a sparse
# checkout of directories.  For example:
#
#   $ cd Trilinos
#   $ ./commonTools/git/sparse_checkout.sh doc packages
#   $ ls
#
# should return just 'commonTools/git/sparse_checkout', 'doc' and
# 'packages'.
#
# After a sparse checkout, commands like 'git status' will run as fast
# as the size of the remaining checkout.
#
# A sparse checkout can be undone with:
#
#   $ git config core.sparsecheckout false
#   $ git read-tree -m -u HEAD
#
# NOTES:
#
# a) Make sure you run this script before you change any files or the
# 'git read-tree ...' command will blow them away (a lesson I learned
# the hard way).  In general, commit your work *before* you run this
# script.
#
# b) Git will keep directories not listed if there are symbolic links
# needed to support files in directories that are listed.

DIRS_FILES_LIST=$@

# Tell Git to allow sparse checkout
git config core.sparsecheckout true

# Set up the list of dirs/files for sparse checkout to keep
SC_FILE=.git/info/sparse-checkout
echo commonTools/git/sparse_checkout.sh > $SC_FILE
for dir_or_file in $DIRS_FILES_LIST ; do
  echo "File/directory to keep: $dir_or_file";
  echo $dir_or_file >> $SC_FILE
done

# Tell git to process the list of dirs/files
git read-tree -m -u HEAD
