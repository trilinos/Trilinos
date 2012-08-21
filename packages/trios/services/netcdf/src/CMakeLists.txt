INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)
INCLUDE(TribitsAddTest)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TriosProcessXDR)


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

# The rpcgen files should be treated as system files to avoid annoying warnings
TriosProcessXDR(${CMAKE_CURRENT_SOURCE_DIR}/netcdf_args.x)


APPEND_SET(NOINSTHEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/netcdf_args.h
  aggregation.h
  create_subchunks.h
  io_timer.h
  NcAttInfo.h
  NcDimInfo.h
  NcFileInfo.h
  NcGroupInfo.h
  NcVarInfo.h
  netcdf_client_private.h
  netcdf_debug.h
)

APPEND_SET(LIB_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/netcdf_args.c
    netcdf_debug.c
    NcFileInfo.cpp
    NcAttInfo.cpp
    NcVarInfo.cpp
    NcDimInfo.cpp
    NcGroupInfo.cpp
    netcdf_client.cpp
    ncmpi_client.cpp
    create_subchunks.cpp
)

############# libtrios_netcdf_client ##########################

TRIBITS_ADD_LIBRARY(
   trios_netcdf_client
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${LIB_SOURCES}
   TESTONLY
)

# trios_netcdf_client is created as a TESTONLY library which means it won't
# get installed.  forced the install here.
INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/libtrios_netcdf_client.a
  DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
)
GLOBAL_SET(${PACKAGE_NAME}_NETCDF_LIBRARIES "trios_netcdf_client")

######################  PROGRAMS  ##################################################

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

APPEND_SET(SVC_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/netcdf_args.c
    netcdf_server.cpp
    netcdf_debug.c
    aggregation.cpp
)

TRIBITS_ADD_EXECUTABLE(
  netcdf-server
  SOURCES ${SVC_SOURCES}
  DEPLIBS ${Trios_LIBRARIES}
  NOEXEPREFIX
)
