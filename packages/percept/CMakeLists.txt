#################
# Percept Package
#################

TRIBITS_PACKAGE(Percept)

#################
# Options
#################

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  PERCEPT_DEBUG
  "Build with debug code enabled."
  ${Trilinos_ENABLE_DEBUG})

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_GEOMETRY
  STK_PERCEPT_HAS_GEOMETRY
  "Build with geometry support."
  OFF)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_LITE
  STK_PERCEPT_LITE
  "Build with geometry support."
  OFF)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ACIS
  HAVE_ACIS
  "Turn off acis support."
  OFF)

#PROPERTY COMPILE_DEFINITIONS STK_PERCEPT_LITE=0 STK_PERCEPT_HAS_GEOMETRY STK_PERCEPT_USE_INTREPID
# change to true trilinos dependency???
#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Intrepid
#  STK_PERCEPT_USE_INTREPID
#  "Build with intrepid support."
#  OFF)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#################
# Add library
#################

#ADD_SUBDIRECTORY(src)
#TRIBITS_ADD_TEST_DIRECTORIES(test)
#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

FILE (GLOB ADAPT_MAIN
     src/adapt/main/*.cpp
)

FILE (GLOB ADAPT_HEADER
     src/adapt/main/*.hpp
)

FILE (GLOB SOURCE

     src/adapt/*.cpp
     src/adapt/markers/*.cpp
     src/adapt/sierra_element/*.cpp

     src/percept/*.cpp
     src/percept/element/intrepid/*.cpp

     src/percept/fixtures/*.cpp
     src/percept/function/*.cpp
     src/percept/function/internal/*.cpp

     src/percept/mesh/gen/*.cpp

     src/percept/mesh/geometry/kernel/*.cpp
     src/percept/mesh/geometry/kernel/search/*.cpp
     src/percept/mesh/geometry/kernel/xfer/*.cpp

     src/percept/mesh/geometry/recovery/*.cpp
     src/percept/mesh/geometry/stk_geom/*.cpp
     src/percept/mesh/geometry/stk_geom/3D/*.cpp
     src/percept/mesh/geometry/volume/*.cpp
     src/percept/mesh/geometry/volume/sierra_only/*.cpp

     src/percept/mesh/mod/smoother/*.cpp

     src/percept/norm/*.cpp

     src/percept/stk_rebalance/*.cpp
     src/percept/stk_rebalance_utils/*.cpp
     src/percept/structured/*.cpp

     src/percept/util/*.cpp
     src/percept/xfer/*.cpp
)

FILE (GLOB HEADER

     src/adapt/*.hpp
     src/adapt/sierra_element/*.hpp
     src/adapt/markers/*.hpp

     src/percept/*.hpp

     src/percept/element/intrepid/*.hpp

     src/percept/fixtures/*.hpp
     src/percept/function/*.hpp
     src/percept/function/internal/*.hpp

     src/percept/mesh/gen/*.hpp
     src/percept/mesh/geometry/kernel/*.hpp
     src/percept/mesh/geometry/kernel/xfer/*.hpp
     src/percept/mesh/geometry/recovery/*.hpp
     src/percept/mesh/geometry/stk_geom/*.hpp
     src/percept/mesh/geometry/stk_geom/3D/*.hpp
     src/percept/mesh/geometry/volume/*.hpp
     src/percept/mesh/geometry/volume/sierra_only/*.hpp
     src/percept/mesh/mod/smoother/*.hpp

     src/percept/norm/*.hpp

     src/percept/stk_rebalance/*.hpp
     src/percept/stk_rebalance_utils/*.hpp
     src/percept/structured/*.hpp

     src/percept/util/*.hpp
     src/percept/xfer/*.hpp
)

TRIBITS_ADD_LIBRARY(
  percept
  HEADERS ${HEADER}
  SOURCES ${SOURCE}
  )

TRIBITS_ADD_EXECUTABLE(
  perceptX
  NOEXEPREFIX
  NOEXESUFFIX
  SOURCES ${ADAPT_MAIN} ${ADAPT_HEADER}
  COMM serial mpi
  )

# Temporary hack until we can figure out how to pull in the defines
# from the config file.
SET_PROPERTY(
  TARGET percept perceptX
  PROPERTY COMPILE_DEFINITIONS STK_PERCEPT_LITE=0 STK_PERCEPT_USE_INTREPID
)

install_executable(perceptX)

TRIBITS_PACKAGE_POSTPROCESS()
