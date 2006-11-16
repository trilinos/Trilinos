#!/bin/sh
#
# Call as
#
#   replace-typeid-name-with-typeName-20061112.sh filename
# 
# and this will replace the construct:
#
#    typeid(anything).name()
#
# with
#
#    typeName(anything)
#
# You can run this on a bunch of source files like:
#
#    find . -name "*pp" -exec \
#      $TRILINOS_HOME/packages/teuchos/refactoring/replace-typeid-name-with-typeName-20061112.sh {} \;
#
# This will cause your code to use the new Teuchos::typeName(...) templated function that
# will use name demangling to print your object types.  Note that this script will also
# replace the versions of typeid(TypeName).name() and will not compile afterward.  For these,
# just manually change it to Teuchos::TypeNameTraits<TypeName>(TypeName)::name().
#

sed --in-place 's/typeid(\(.\+\)).name()/typeName(\1)/g' $1
