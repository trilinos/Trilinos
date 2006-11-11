#!/bin/sh
#
# Call as
#
#   replace-OSTab-operator-20061110.sh filename
# 
# and this will replace the construct:
#
#    OSTab(anything)()
#
# with
#
#    OSTab(anything).o()
#
# You can run this on a bunch of source files like:
#
#    find . -name "*pp" -exec \
#      $TRILINOS_HOME/packages/teuchos/refactoring/replace-OSTab-operator-20061110.sh {} \;
#
# This is to fix a bug with the gcc 3.3.3 compiler that Ken
# Stanly is trying to use that is on a system at the company
# StarP.
#

sed --in-place 's/OSTab(\(.\+\))()/OSTab(\1)\.o()/g' $1
