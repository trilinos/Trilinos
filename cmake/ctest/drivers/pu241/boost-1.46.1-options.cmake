#
# Set up TPL pointers for boost 1.46.1
#

SET(TPL_ENABLE_Boost ON CACHE BOOL "")
SET(Boost_INCLUDE_DIRS /opt/gcc-4.5.1/tpls/boost-1.46.1/include CACHE FILEPATH "")
SET(BoostLib_LIBRARY_NAMES "boost_thread;boost_program_options" CACHE STRING  "")
SET(BoostLib_INCLUDE_DIRS ${Boost_INCLUDE_DIRS}  CACHE FILEPATH "")
SET(BoostLib_LIBRARY_DIRS /opt/gcc-4.5.1/tpls/boost-1.46.1/lib  CACHE FILEPATH "")

# Must force this variable in the Dakota configure process
SET(BOOST_INCLUDEDIR ${Boost_INCLUDE_DIRS}  CACHE PATH  "" FORCE)
