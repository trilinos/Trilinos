
#
# A) Define the package
#

TRIBITS_PACKAGE_DECL( SEACAS ENABLE_SHADOWING_WARNINGS )

TRIBITS_ADD_DEBUG_OPTION()
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(
  SEACASExodus_ENABLE_THREADSAFE
  EXODUS_THREADSAFE
  "Enable thread-safe version of Exodus."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  SEACASIoss_ENABLE_THREADSAFE
  IOSS_THREADSAFE
  "Enable thread-safe version of Ioss."
  OFF
  )

function(disable_if_no_fortran package)
IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
  IF (${PROJECT_NAME}_ENABLE_SEACAS${package})
    MESSAGE("-- "
      "WARNING: Setting ${PROJECT_NAME}_ENABLE_SEACAS${package}=OFF"
      " even though it was set to ${${PROJECT_NAME}_ENABLE_SEACAS${package}}"
      " because ${PROJECT_NAME}_ENABLE_Fortran=OFF!"
      )
  ELSEIF("${${PROJECT_NAME}_ENABLE_SEACAS${package}}" STREQUAL "")
     MESSAGE("-- "
      "NOTE: Setting ${PROJECT_NAME}_ENABLE_SEACAS${package}=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran=OFF!"
       )
  ENDIF()
  SET(${PROJECT_NAME}_ENABLE_SEACAS${package} OFF PARENT_SCOPE)
ENDIF()
endfunction()

disable_if_no_fortran(Mapvarlib)
disable_if_no_fortran(Exodus_for)
disable_if_no_fortran(ExoIIv2for32)
disable_if_no_fortran(Supes)
disable_if_no_fortran(Suplib)
disable_if_no_fortran(PLT)
disable_if_no_fortran(Blot)
disable_if_no_fortran(Fastq)
disable_if_no_fortran(SVDI)
disable_if_no_fortran(Algebra)
disable_if_no_fortran(Exotxt)
disable_if_no_fortran(Gjoin)
disable_if_no_fortran(Gen3D)
disable_if_no_fortran(Genshell)
disable_if_no_fortran(Grepos)
disable_if_no_fortran(Explore)
disable_if_no_fortran(Mapvar)
disable_if_no_fortran(Mapvar-kd)
disable_if_no_fortran(Numbers)
disable_if_no_fortran(Txtexo)
disable_if_no_fortran(Ex2ex1v2)
disable_if_no_fortran(Ex1ex2v2)

IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
# For some reason, this variable is undefined if fortran is disabled and it causes a cmake error.
# define it to a nonsense variable to avoid error.
  IF (NOT CMAKE_Fortran_LINK_EXECUTABLE)
    MESSAGE("-- "
      "WARNING: Setting CMAKE_Fortran_LINK_EXECUTABLE to a random value to avoid CMake error"
      )
    SET(CMAKE_Fortran_LINK_EXECUTABLE empty)
  ENDIF()
ENDIF()

ASSERT_DEFINED(TPL_ENABLE_Pthread)
IF (SEACASExodus_ENABLE_THREADSAFE AND NOT TPL_ENABLE_Pthread)
  MESSAGE(FATAL_ERROR "You set SEACASExodus_ENABLE_THREADSAFE=ON, but TriBITs support for Pthread(s) is not enabled (TPL_ENABLE_Pthread=OFF).  This is not allowed.  Please enable Pthreads before attempting to enable the thread-safe Exodus support.")
ENDIF ()
IF (SEACASIoss_ENABLE_THREADSAFE AND NOT TPL_ENABLE_Pthread)
  MESSAGE(FATAL_ERROR "You set SEACASIoss_ENABLE_THREADSAFE=ON, but TriBITs support for Pthread(s) is not enabled (TPL_ENABLE_Pthread=OFF).  This is not allowed.  Please enable Pthreads before attempting to enable the thread-safe IOSS support.")
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

