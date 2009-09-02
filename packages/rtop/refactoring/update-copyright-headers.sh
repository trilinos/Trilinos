#!/usr/bin/env sh
#
# This script automatically updates the copyright headers for
# for all files.
#
# To run this script, you must have the variable $TRILINOS_HOME set
# to the base Trilinos directory so that
#
#    $TRILINOS_HOME/packages/rtop/refactoring
#
# is where this file is.  Trilinos/commonTools/refactoring must also
# be added to your path.
#
# This script can be run for the base directory $TRILINOS_HOME/package/rtop
# over and over again without any ill side effects.
#

find . -name "*.h" -exec $TRILINOS_HOME/commonTools/refactoring/update-copyright-header.pl \
  $TRILINOS_HOME/packages/rtop/refactoring/c++-copyright-header.txt {} {} \;

find . -name "*.hpp" -exec $TRILINOS_HOME/commonTools/refactoring/update-copyright-header.pl \
  $TRILINOS_HOME/packages/rtop/refactoring/c++-copyright-header.txt {} {} \;

find . -name "*.cpp" -exec $TRILINOS_HOME/commonTools/refactoring/update-copyright-header.pl \
  $TRILINOS_HOME/packages/rtop/refactoring/c++-copyright-header.txt {} {} \;

find . -name "Makefile.am" -exec $TRILINOS_HOME/commonTools/refactoring/update-copyright-header.pl \
  $TRILINOS_HOME/packages/rtop/refactoring/scirpt-copyright-header.txt {} {} \;
