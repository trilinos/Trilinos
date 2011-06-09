INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( SEACAS ENABLE_SHADOWING_WARNINGS )

PACKAGE_ADD_DEBUG_OPTION()
PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

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
IF (NOT TPL_ENABLE_MPI)
INSTALL(TARGETS ${target} COMPONENT ${PACKAGE_NAME} RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR})
ENDIF()
endfunction()

SET(CMAKE_EXECUTABLE_SUFFIX "")
ADD_SUBDIRECTORY(libraries)

IF (${PACKAGE_NAME}_ENABLE_APPLICATIONS)
ADD_SUBDIRECTORY(applications)
#ADD_SUBDIRECTORY(scripts)
ENDIF()

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()

