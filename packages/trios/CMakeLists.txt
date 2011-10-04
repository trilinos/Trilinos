INCLUDE(PackageMacros)
INCLUDE(PackageLibraryMacros)
INCLUDE(AddOptionAndDefine)

###########################################################
# A) Define the package
#

PACKAGE( Trios ENABLE_SHADOWING_WARNINGS )

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()


###########################################################
# B) Set up package-specific options
#

# Checks for various types of timers
INCLUDE(TriosProbeTimers)
INCLUDE(TriosProbeNetwork)

IF (TPL_ENABLE_Netcdf)
  ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_NETCDF_SERVICE
    TRIOS_ENABLE_NETCDF_SERVICE
    "Build libraries for netcdf service."
    OFF )
ENDIF (TPL_ENABLE_Netcdf)

IF (TPL_ENABLE_Pnetcdf)
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_PNETCDF_SERVICE
   TRIOS_ENABLE_PNETCDF_SERVICE
   "Build libraries for PnetCDF service."
   OFF )
ENDIF (TPL_ENABLE_Pnetcdf)

ADD_OPTION_AND_DEFINE( ${PACKAGE_NAME}_ENABLE_TRACING
  TRIOS_ENABLE_TRACING
  "Enable the tracing API"
  OFF
)

# Check for an RDMA transport (needed for Nessie)
IF (${PACKAGE_NAME}_ENABLE_Portals OR
    ${PACKAGE_NAME}_ENABLE_CrayPortals OR
    ${PACKAGE_NAME}_ENABLE_InfiniBand OR
    ${PACKAGE_NAME}_ENABLE_Gemini OR 
    ${PACKAGE_NAME}_ENABLE_LUC)

  SET(${PACKAGE_NAME}_HAVE_RDMA_TRANSPORT 1)        
ELSE()
  SET(${PACKAGE_NAME}_HAVE_RDMA_TRANSPORT 0)
ENDIF()



# Create the Trios_config.h file in ${CMAKE_BINARY_DIR}
PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

# Copy the ProcessXDR.cmake file
CONFIGURE_FILE(
  ${PACKAGE_SOURCE_DIR}/cmake/ProcessXDR.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/ProcessXDR.cmake
  COPYONLY
)

# Copy the ProcessXDR.cmake file
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


ADD_SUBDIRECTORY(libraries)
ADD_SUBDIRECTORY(programs)
ADD_SUBDIRECTORY(services)

PACKAGE_ADD_TEST_DIRECTORIES(tests)
PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)



###########################################################
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()


# Install cmake scripts so others can use them
INSTALL(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/ProcessXDR.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/UseTrios.cmake
  DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}/cmake/${PACKAGE_NAME}"
)

MESSAGE(STATUS "Trios_TPL_LIBRARIES = ${TPL_LIBRARIES}")
