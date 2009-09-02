
#
# Macro that sets a variable name both in the current scope and the
# parent scope.
#
# It turns out that when you call ADD_SUBDIRECTORY(someDir) that CMake
# actaully creates a copy of all of the regular non-cache varaibles in
# the current scope in order to create a new set of variables for the
# CMakeLists.txt file in 'someDir'.  This means that if you call
# SET(SOMEVAR Blah PARENT_SCOPE) that it will not affect the value of
# SOMEVAR in the current scope.  This macro therefore is designed to
# set the value of the variable in the current scope and the parent
# scope in one shot.
#
# Global variables are different.  When you move to a subordinate
# CMakeLists.txt file, a local copy of the variable is *not* created.
# If you set the value name locally, it will shadow the global
# variable.  However, if you set the globlal value with SET(SOMEVAR
# someValue CACHE INTERNAL ""), then the value will get changed in the
# current subordinate scope and in all parent scopes.
#

MACRO(DUAL_SCOPE_SET VARNAME)
  SET(${VARNAME} ${ARGN} PARENT_SCOPE)
  SET(${VARNAME} ${ARGN})
ENDMACRO()
