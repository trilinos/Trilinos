INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TriosProcessXDR)

# INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

SET(NNTI_HEADERS "")
SET(NNTI_NONINSTALLHEADERS "")
SET(NNTI_SOURCES "")

INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_BINARY_DIR})

############# Special exception for rpcgen extra commas in enum ##########
IF (CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)


############# COMMON ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/..)


TRIBITS_SUBPACKAGE(nnti)


############# NNTI ##########################

# The rpcgen files should be treated as system files to avoid annoying warnings
# INCLUDE_DIRECTORIES(SYSTEM ${CMAKE_CURRENT_BINARY_DIR})

TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/Trios_nnti_xdr.x)

APPEND_SET(NNTI_NOINSTHEADERS
  nnti_gni.h
  nnti_ib.h
  nnti_mpi.h
  nnti_internal.h
  nnti_ptls.h
  nnti_utils.h
  buffer_queue.h
)

APPEND_SET(NNTI_HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nnti_xdr.h
  Trios_nnti_fprint_types.h
  Trios_nnti.h
)

APPEND_SET(NNTI_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nnti_xdr.c
  fprint_types.cpp
  nnti.c
  nnti_utils.c
  buffer_queue.cpp
)

SET(TRIOS_SUPPORTED_NETWORK_FOUND 0)
IF (${PACKAGE_NAME}_ENABLE_Portals OR ${PACKAGE_NAME}_ENABLE_CrayPortals)
  APPEND_SET(NNTI_SOURCES nnti_ptls.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
#IF (${PACKAGE_NAME}_ENABLE_CrayPortals)
#  APPEND_SET(NNTI_SOURCES nnti_ptls.cpp)
#  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
#ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_InfiniBand)
  APPEND_SET(NNTI_SOURCES nnti_ib.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_LUC)
  APPEND_SET(NNTI_SOURCES nnti_luc.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_Gemini)
  APPEND_SET(NNTI_SOURCES nnti_gni.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (${PARENT_PACKAGE_NAME}_ENABLE_MPI)
  APPEND_SET(NNTI_SOURCES nnti_mpi.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (NOT TRIOS_SUPPORTED_NETWORK_FOUND)
   message(FATAL "Did not find a supported network protocol. ")
ENDIF ()

############# libtrios_nnti.a ##########################

TRIBITS_ADD_LIBRARY(
   trios_nnti
   HEADERS ${NNTI_HEADERS}
   NOINSTALLHEADERS ${NNTI_NOINSTALLHEADERS}
   SOURCES ${NNTI_SOURCES}
   DEPLIBS ${DEPLIBS}
)


TRIBITS_SUBPACKAGE_POSTPROCESS()
