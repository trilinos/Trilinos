INCLUDE(TribitsLibraryMacros)
INCLUDE(CheckIncludeFile)

TRIBITS_SUBPACKAGE(Ioss)
TRIBITS_CONFIGURE_FILE(SEACASIoss_config.h)

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/SEACASIoss_config.h)

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/Ioss_*.h)
APPEND_GLOB(SOURCES ${DIR}/*.C)

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_BINARY_DIR}"
  )

TRIBITS_ADD_LIBRARY(
  Ioss
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
)

#
# If built inside Trilinos, then the Ioss library is the subpackage
# and the other libraries are built via the add_subdirectory calls
#
# If built as the standalone SEACAS project, then Ioss is the package
# and the other libraries are subpackages which are listed in
# cmake/Dependencies.cmake

if(${CMAKE_PROJECT_NAME} STREQUAL "SEACAS")
  TRIBITS_PACKAGE_POSTPROCESS()
else()
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_SEACASExodus)
  IF (${${PACKAGE_NAME}_ENABLE_SEACASExodus})
    ADD_SUBDIRECTORY(exodus)
    ADD_SUBDIRECTORY(exo_fac)
    ADD_SUBDIRECTORY(exo_fpp)

    IF (SEACAS_Exodus_ENABLE_MPI)
      ADD_SUBDIRECTORY(exo_par)
    ENDIF()
  ENDIF()

  ASSERT_DEFINED(${CMAKE_PROJECT_NAME}_ENABLE_Pamgen)
  IF (${CMAKE_PROJECT_NAME}_ENABLE_Pamgen)
    ADD_SUBDIRECTORY(pamgen)
  ENDIF()

  ADD_SUBDIRECTORY(generated)
  ADD_SUBDIRECTORY(heartbeat)
  ADD_SUBDIRECTORY(transform)
  ADD_SUBDIRECTORY(init)
  ADD_SUBDIRECTORY(main)
  TRIBITS_ADD_TEST_DIRECTORIES(utest)
  TRIBITS_SUBPACKAGE_POSTPROCESS()
endif()


