#!/bin/sh
#
# Call as
#
#   add-semicolon-after-composition-members-macros-20061116.sh filename
# 
# and this will replace the construct:
#
#    STANDARD_.*COMPOSITION_MEMBERS(anything)
#
# with
#
#    STANDARD_.*COMPOSITION_MEMBERS(anything);
#
# You can run this on a bunch of source files like:
#
#    find . -name "*pp" -exec \
#      $TRILINOS_HOME/packages/teuchos/refactoring/add-semicolon-after-composition-members-macros-20061116.sh {} \;
#
# This is needed to update code after I refactored these macros
# to require a semicolon (see Bug 2897).
#

sed --in-place 's/\(STANDARD_.*COMPOSITION_MEMBERS\)(\(.\+\))/\1(\2);/g' $1
