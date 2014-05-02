#
# Set up TPL pointers for boost 1.49.0
#

SET(BOOST_ROOT /opt/intel-12.0.8/tpls/boost-1.49.0)
SET(TPL_ENABLE_Boost ON CACHE BOOL "")
SET(Boost_INCLUDE_DIRS ${BOOST_ROOT}/include CACHE FILEPATH "")
SET(BoostLib_LIBRARY_NAMES "boost_thread;boost_program_options" CACHE STRING  "")
SET(BoostLib_INCLUDE_DIRS ${Boost_INCLUDE_DIRS}  CACHE FILEPATH "")
SET(BoostLib_LIBRARY_DIRS ${BOOST_ROOT}/lib  CACHE FILEPATH "")

# Must force this variable in the Dakota configure process
SET(BOOST_INCLUDEDIR ${Boost_INCLUDE_DIRS}  CACHE PATH  "" FORCE)
