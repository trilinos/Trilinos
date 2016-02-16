INCLUDE(TribitsLibraryMacros)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/exodusII_cfg.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/exodusII_cfg.h
  @ONLY IMMEDIATE)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_BINARY_DIR}/include"
)

FILE(GLOB SOURCES src/*.c)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

IF (TPL_Netcdf_Enables_Netcdf4)
   ADD_DEFINITIONS(-DENABLE_NETCDF4)
ENDIF()

SET(HEADERS
  include/exodusII.h
  include/exodusII_par.h
  ${CMAKE_CURRENT_BINARY_DIR}/include/exodusII_meta.h
)

#####
# Create 'exodusII_meta.h' include file.
#####

# A function used to create 1/0 definitions.
# If a variable is set, then '1' is returned. Otherwise, '0' is
# returned. The ret_val is prepended with 'EX'
# when feature is true, which is used to generate exodusII_meta.h.
FUNCTION(is_enabled feature ret_val)
  IF(${feature})
    SET("EX_${ret_val}" 1 PARENT_SCOPE)
  ELSE()
    SET("EX_${ret_val}" 0 PARENT_SCOPE)
  ENDIF(${feature})
ENDFUNCTION()

is_enabled(TPL_Netcdf_Enables_Netcdf4 HAS_NETCDF4)
is_enabled(TPL_Netcdf_Enables_PNetcdf HAS_PNETCDF)
is_enabled(TPL_Netcdf_PARALLEL HAS_PARALLEL)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/exodusII_meta.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/exodusII_meta.h @ONLY)

IF (NOT TPL_Netcdf_PARALLEL)
  LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/ex_create_par.c")
  LIST(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/ex_open_par.c")
  LIST(REMOVE_ITEM HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/exodusII_par.h")
ENDIF()

TRIBITS_ADD_LIBRARY(
	exodus
	HEADERS	${HEADERS}
	SOURCES ${SOURCES}
)

TRIBITS_ADD_TEST_DIRECTORIES(test)

