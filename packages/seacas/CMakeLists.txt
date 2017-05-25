
#
# A) Define the package
#

TRIBITS_PACKAGE_DECL( SEACAS ENABLE_SHADOWING_WARNINGS )
TRIBITS_PACKAGE_DECL( SEACAS ENABLE_STRONG_C_COMPILE_WARNINGS )
TRIBITS_PACKAGE_DECL( SEACAS ENABLE_STRONG_CXX_COMPILE_WARNINGS )

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(
  SEACASExodus_ENABLE_THREADSAFE
  EXODUS_THREADSAFE
  "Enable thread-safe version of Exodus."
  OFF
  )

ASSERT_DEFINED(TPL_ENABLE_Pthread)
IF (EXODUS_ENABLE_THREADSAFE AND NOT TPL_ENABLE_Pthread)
  MESSAGE(FATAL_ERROR "You set EXODUS_ENABLE_THREADSAFE=ON, but TriBITs support for Pthread(s) is not enabled (TPL_ENABLE_Pthread=OFF).  This is not allowed.  Please enable Pthreads before attempting to enable the thread-safe Exodus support.")
ENDIF ()

# C) Add the libraries, tests, and examples
#- Add executable that will be installed in bin dir...
#
#  ADD_BIN_EXECUTABLE(<executable>)
#
function(install_executable target)
  INSTALL(TARGETS ${target} COMPONENT ${PACKAGE_NAME} RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR})
endfunction()

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

ADD_SUBDIRECTORY(scripts)

TRIBITS_EXCLUDE_FILES(
  doc
  )
  
#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

