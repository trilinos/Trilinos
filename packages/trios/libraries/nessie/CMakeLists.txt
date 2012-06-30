
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TriosProcessXDR)

############# Process Subdirectories ###################
#ADD_SUBDIRECTORY(common)


# INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_BINARY_DIR})

############# Special exception for rpcgen extra commas in enum ##########
IF (CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)

############# COMMON ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/common)

# The rpcgen files should be treated as system files to avoid annoying warnings
# INCLUDE_DIRECTORIES(SYSTEM ${CMAKE_CURRENT_BINARY_DIR})

TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/common/Trios_nnti_xdr.x)
TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/common/Trios_nssi_types_xdr.x)
TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/common/nssi_service_args.x)


APPEND_SET(NOINSTHEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/nssi_service_args.h
  common/nnti.h
  common/nnti_gni.h
  common/nnti_ib.h
  common/nnti_mpi.h
  common/nnti_internal.h
  common/nnti_ptls.h
  common/nnti_utils.h
  common/nssi_opcodes.h
  common/nssi_trace.h
  common/buffer_queue.h
)

APPEND_SET(HEADERS
  common/Trios_nssi_rpc.h
  common/Trios_nssi_types.h
  common/Trios_nssi_debug.h
  common/Trios_nssi_fprint_types.h
  common/Trios_nssi_xdr.h
  common/Trios_nssi_TypeManip.h
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nnti_xdr.h
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nssi_types_xdr.h
)


APPEND_SET(SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nnti_xdr.c
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nssi_types_xdr.c
  ${CMAKE_CURRENT_BINARY_DIR}/nssi_service_args.c
  common/fprint_types.cpp
  common/nnti.c
  common/nnti_utils.c
  common/nssi_common.cpp
  common/nssi_debug.c
  common/nssi_xdr.cpp
  common/buffer_queue.cpp
)

SET(TRIOS_SUPPORTED_NETWORK_FOUND 0)
IF (${PACKAGE_NAME}_ENABLE_Portals OR ${PACKAGE_NAME}_ENABLE_CrayPortals)
  APPEND_SET(SOURCES common/nnti_ptls.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
#IF (${PACKAGE_NAME}_ENABLE_CrayPortals)
#  APPEND_SET(SOURCES common/nnti_ptls.cpp)
#  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
#ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_InfiniBand)
  APPEND_SET(SOURCES common/nnti_ib.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_LUC)
  APPEND_SET(SOURCES common/nnti_luc.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_Gemini)
  APPEND_SET(SOURCES common/nnti_gni.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (${PACKAGE_NAME}_ENABLE_MPI)
  APPEND_SET(SOURCES common/nnti_mpi.cpp)
  SET(TRIOS_SUPPORTED_NETWORK_FOUND 1)
ENDIF ()
IF (NOT TRIOS_SUPPORTED_NETWORK_FOUND)
   message(FATAL "Did not find a supported network protocol. ")
ENDIF ()


############# CLIENT ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/client)

APPEND_SET(HEADERS
  client/Trios_nssi_client.h
  client/Trios_nssi_request.h
)

APPEND_SET(SOURCES
  client/nssi_client.cpp
)


############# SERVER ##########################


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/server)

APPEND_SET(HEADERS
  server/Trios_nssi_server.h
)


APPEND_SET(SOURCES
  server/nssi_server.cpp
)


############# libnssi ##########################

TRIBITS_ADD_LIBRARY(
   trios_nssi
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${SOURCES}
   DEPLIBS ${DEPLIBS} trios_support
)

