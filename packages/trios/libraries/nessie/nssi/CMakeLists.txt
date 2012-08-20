INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TriosProcessXDR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

SET(NSSI_HEADERS "")
SET(NSSI_NONINSTALLHEADERS "")
SET(NSSI_SOURCES "")

INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_SOURCE_DIR})

############# Special exception for rpcgen extra commas in enum ##########
IF (CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
ENDIF (CMAKE_COMPILER_IS_GNUCXX)


############# COMMON ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../nnti)


TRIBITS_SUBPACKAGE(nssi)


############# NSSI ##########################

# The rpcgen files should be treated as system files to avoid annoying warnings
# INCLUDE_DIRECTORIES(SYSTEM ${CMAKE_CURRENT_BINARY_DIR})

TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/Trios_nssi_types_xdr.x)
TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/nssi_service_args.x)

APPEND_SET(NSSI_NOINSTHEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/nssi_service_args.h
  nssi_opcodes.h
  nssi_trace.h
)

APPEND_SET(NSSI_HEADERS
  Trios_nssi_rpc.h
  Trios_nssi_types.h
  Trios_nssi_debug.h
  Trios_nssi_fprint_types.h
  Trios_nssi_xdr.h
  Trios_nssi_TypeManip.h
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nssi_types_xdr.h
)

APPEND_SET(NSSI_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/Trios_nssi_types_xdr.c
  ${CMAKE_CURRENT_BINARY_DIR}/nssi_service_args.c
  fprint_types.cpp
  nssi_common.cpp
  nssi_debug.c
  nssi_xdr.cpp
)

############# CLIENT ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../client)

APPEND_SET(NSSI_HEADERS
  ../client/Trios_nssi_client.h
  ../client/Trios_nssi_request.h
)

APPEND_SET(NSSI_SOURCES
  ../client/nssi_client.cpp
)

############# SERVER ##########################

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../server)

APPEND_SET(NSSI_HEADERS
  ../server/Trios_nssi_server.h
)

APPEND_SET(NSSI_SOURCES
  ../server/nssi_server.cpp
)

############# libtrios_nssi.a ##########################

TRIBITS_ADD_LIBRARY(
   trios_nssi
   HEADERS ${NSSI_HEADERS}
   NOINSTALLHEADERS ${NSSI_NOINSTALLHEADERS}
   SOURCES ${NSSI_SOURCES}
   DEPLIBS ${DEPLIBS}
)


TRIBITS_SUBPACKAGE_POSTPROCESS()
