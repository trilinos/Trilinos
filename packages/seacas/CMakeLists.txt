
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

set(${PROJECT_NAME}_C_Standard_DEFAULT c99)

SET(${PROJECT_NAME}_CMAKE_CXX_STANDARD_DEFAULT 17)
SET(CMAKE_CXX_STANDARD 17)
SET(${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED "(14|17|20)")
ADVANCED_SET(CMAKE_CXX_STANDARD ${${PROJECT_NAME}_CMAKE_CXX_STANDARD_DEFAULT}
  CACHE STRING
  "C++ standard number with values ${${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED} (default ${${PROJECT_NAME}_CMAKE_CXX_STANDARD_DEFAULT})")
IF (NOT CMAKE_CXX_STANDARD MATCHES "^${${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED}$")
  MESSAGE(FATAL_ERROR
    "CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} is not in the allowed set "
    "${${PROJECT_NAME}_CMAKE_CXX_STANDARDS_ALLOWED}")
ENDIF ()

if ( CYGWIN )
  # This seems to be needed to get `realpath` defined in the IOSS build...
  set(CMAKE_CXX_EXTENSIONS ON)
endif()

set(SEACASExodus_ENABLE_STATIC_DEFAULT "OFF")
set(SEACASExodus_ENABLE_SHARED_DEFAULT "OFF")

# On windows (and probably other systems), we only want a single
# exodus library and not both static and shared.
if (NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  if (${CMAKE_PROJECT_NAME} STREQUAL "Seacas")
    if (BUILD_SHARED_LIBS)
      set(SEACASExodus_ENABLE_STATIC_DEFAULT "ON")
    else()
      set(SEACASExodus_ENABLE_SHARED_DEFAULT "ON")
    endif()
  endif()

  TRIBITS_ADD_OPTION_AND_DEFINE(
    SEACASExodus_ENABLE_STATIC
    EXODUS_STATIC
    "Build static version of Exodus in addition to shared version."
    ${SEACASExodus_ENABLE_STATIC_DEFAULT}
    )
  TRIBITS_ADD_OPTION_AND_DEFINE(
    SEACASExodus_ENABLE_SHARED
    EXODUS_SHARED
    "Build shared version of Exodus in addition to static version."
    ${SEACASExodus_ENABLE_SHARED_DEFAULT}
    )
endif()

TRIBITS_ADD_OPTION_AND_DEFINE(
  SEACASIoss_ENABLE_THREADSAFE
  IOSS_THREADSAFE
  "Enable thread-safe version of Ioss."
  OFF
  )

ADVANCED_SET(${CMAKE_PROJECT_NAME}_ENABLE_SEACASIoss_Applications ON CACHE BOOL
								  "Build the io_shell, io_info, io_modify and other ioss executables")

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
  IF (TARGET ${target})
    INSTALL(TARGETS ${target} COMPONENT ${PACKAGE_NAME} RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_RUNTIME_DIR})
  ENDIF()
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
