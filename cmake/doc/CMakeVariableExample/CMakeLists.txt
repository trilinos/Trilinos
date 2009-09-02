# This simple project demonstrates how CMake handles variables.  Just
# cd to a different scratch directory and type 'cmake --debug-output
# THISDIR' and you will see what happens.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(CMakeVaraibleExample NONE)

# Create a regular local variable that exists in the current (base)
# scope
SET(SOMEVAR baseName)
MESSAGE("Value of SOMEVAR in parent is '${SOMEVAR}'")

# Create a special global variable that will span all scopes
SET(SOMEGLOBALVAR baseName CACHE INTERNAL "")
MESSAGE("Value of SOMEGLOBALVAR in parent is '${SOMEGLOBALVAR}'")

# Move to a subordinate CMakeLists.txt file.  This will create a new
# scope.  In the process, all of the regular local varaibles like
# SOMEVAR will be copied and placed into a new scope.  However, the
# global variables like SOMEGLOBALVAR will not be copied but are still
# accessible in the subordinate scope.

ADD_SUBDIRECTORY(subdir)

# This will print the value of SOMEVAR set above into the parent scope
# (i.e. this scope) and not what was set in the local copy in the
# subdirs/CMakeLists.txt scope.
MESSAGE("Value of SOMEVAR in parent is '${SOMEVAR}'")

# This will print the value of SOMEGLOBALVAR that was set from
# subdir/CMakeLists.txt which should be 'newName3'.
MESSAGE("Value of SOMEGLOBALVAR in parent is '${SOMEGLOBALVAR}'")
