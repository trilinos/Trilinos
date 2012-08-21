INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddOptionAndDefine)

###########################################################
# A) Define the package
#

TRIBITS_PACKAGE_DECL( Trios ENABLE_SHADOWING_WARNINGS )

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()


###########################################################
# B) Set up package-specific options
#

set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

# Various Trios checks
INCLUDE(TriosProbeTimers)
# Checks for a supported interconnect
INCLUDE(TriosProbeNetwork)
# Checks for a POSIX semaphore implementation
INCLUDE(TriosProbeSemaphores)
# Check XDR implementation
INCLUDE(TriosProbeXDR)

IF (TPL_ENABLE_Netcdf)
  TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_NETCDF_SERVICE
    TRIOS_ENABLE_NETCDF_SERVICE
    "Build libraries for netcdf service."
    OFF )
ENDIF (TPL_ENABLE_Netcdf)

IF (TPL_ENABLE_Pnetcdf)
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PNETCDF_SERVICE
   TRIOS_ENABLE_PNETCDF_SERVICE
   "Build libraries for PnetCDF service."
   OFF )
ENDIF (TPL_ENABLE_Pnetcdf)

TRIBITS_ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_TRACING
  TRIOS_ENABLE_TRACING
  "Enable the tracing API"
  OFF
)

# Check for an RDMA transport (needed for Nessie)
IF (${PACKAGE_NAME}_ENABLE_Portals OR
    ${PACKAGE_NAME}_ENABLE_CrayPortals OR
    ${PACKAGE_NAME}_ENABLE_InfiniBand OR
    ${PACKAGE_NAME}_ENABLE_Gemini OR
    ${PACKAGE_NAME}_ENABLE_MPI OR
    ${PACKAGE_NAME}_ENABLE_LUC)

  SET(${PACKAGE_NAME}_HAVE_RDMA_TRANSPORT 1)
ELSE()
  SET(${PACKAGE_NAME}_HAVE_RDMA_TRANSPORT 0)
ENDIF()



# Create the Trios_config.h file in ${CMAKE_BINARY_DIR}
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

# Copy the TriosProcessXDR.cmake file
CONFIGURE_FILE(
  ${PACKAGE_SOURCE_DIR}/cmake/TriosProcessXDR.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/TriosProcessXDR.cmake
  COPYONLY
)

# Copy the UseTrios.cmake file
CONFIGURE_FILE(
  ${PACKAGE_SOURCE_DIR}/cmake/UseTrios.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/UseTrios.cmake
  COPYONLY
)


###########################################################
# C) Add the libraries, tests, and examples
#

# Need this include so codes can find the config file
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()


###########################################################
# D) Do standard postprocessing
#

# Install cmake scripts so others can use them
INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/TriosProcessXDR.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/UseTrios.cmake
  DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
)

MESSAGE(STATUS "Trios_TPL_LIBRARIES = ${TPL_LIBRARIES}")


TRIBITS_PACKAGE_POSTPROCESS()

#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

