INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE_DECL( SEACAS ENABLE_SHADOWING_WARNINGS )
TRIBITS_PACKAGE_DECL( SEACAS ENABLE_STRONG_C_COMPILE_WARNINGS )
TRIBITS_PACKAGE_DECL( SEACAS ENABLE_STRONG_CXX_COMPILE_WARNINGS )

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#
 if (DEFINED TPL_Netcdf_PARALLEL)
 else()
   ASSERT_DEFINED(SEACASExodus_ENABLE_MPI)
   SET(TPL_Netcdf_PARALLEL ${SEACASExodus_ENABLE_MPI})
 endif()

#
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

