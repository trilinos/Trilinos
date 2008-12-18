
#
# Set platform-specific options needed by the platform independent
# code
#

SET (CTEST_CVS_COMMAND "cvs -q -z3")
SET(CTEST_CMAKE_COMMAND /usr/local/bin/cmake)

#
# Read in the platform-independent options
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestSupport.cmake")


#
# Update or add to platform-specific options
#


SET(CTEST_INITIAL_CACHE
"
${CTEST_INITIAL_CACHE}
MEMORYCHECK_COMMAND:FILEPATH=/usr/bin/valgrind
MAKECOMMAND:STRING=make -j8 -i
"
)
