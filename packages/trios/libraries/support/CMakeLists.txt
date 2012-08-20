
INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)

TRIBITS_SUBPACKAGE(support)

# A) Package-specific configuration options
#

# Initializations
SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")
SET(MODULES "")


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

#INCLUDE_DIRECTORIES(BEFORE SYSTEM ${${PARENT_PACKAGE}_BINARY_DIR})

#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/..)

message(STATUS "${${PARENT_PACKAGE_NAME}_BINARY_DIR}/Trios_config.h")

APPEND_SET(HEADERS
  ${${PARENT_PACKAGE_NAME}_BINARY_DIR}/Trios_config.h
)

############# LOGGER ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/logger)

APPEND_SET(HEADERS
  logger/Trios_logger.h
)

APPEND_SET(SOURCES
  logger/logger.cpp
)


############# THREADS ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/threads)

APPEND_SET(HEADERS
  threads/Trios_threads.h
  threads/Trios_threads_types.h
)

APPEND_SET(SOURCES
  threads/trios_threads.cpp
)

############# TIMER ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/timer)

APPEND_SET(HEADERS
  timer/Trios_timer.h
)

APPEND_SET(SOURCES
  timer/timer.c
  timer/Trios_Timer.cpp
)


############# SIGNAL ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/signal)

APPEND_SET(HEADERS
  signal/Trios_signal.h
)

APPEND_SET(SOURCES
  signal/signal.cpp
)


########### Trace  ###############



INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/trace)

APPEND_SET(HEADERS
  trace/Trios_trace.h
)

APPEND_SET(NOINSTALLHEADERS
  trace/SDDF.h
  trace/TraceFile.h
  trace/TraceGroup.h
)

APPEND_SET(SOURCES
  trace/trace.cc
)

# See if PABLO is enabled
IF (${PACKAGE_NAME}_ENABLE_Pablo)
  APPEND_SET(SOURCES trace/SDDF.cc)
ENDIF (${PACKAGE_NAME}_ENABLE_Pablo)



########################### TRIOS SUPPORT LIBRARY ###############

TRIBITS_ADD_LIBRARY(
   trios_support
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${SOURCES}
   DEPLIBS ${DEPLIBS}
   )


########### install targets ###############

# INSTALL(FILES ${MODULES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

TRIBITS_SUBPACKAGE_POSTPROCESS()
