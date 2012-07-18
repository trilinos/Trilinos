INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE_DECL( SEACAS ENABLE_SHADOWING_WARNINGS )

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# B) Set up package-specific options
#

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

IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Nemspread AND ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}Nemslice)
  ADD_SUBDIRECTORY(scripts)
ENDIF()

TRIBITS_EXCLUDE_FILES(
  doc
  )
  
#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

