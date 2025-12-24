TRIBITS_PACKAGE(magistrate)
INCLUDE(FetchContent)

OPTION(${PACKAGE_NAME}_DOWNLOAD_SOURCE "If \"ON\" then magistrate source code will be downloaded from the internet." OFF)
IF(${PROJECT_NAME}_ENABLE_magistrate)
  IF(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/magistrate")
    MESSAGE(STATUS "Enabling inline source build of magistrate = ${HAVE_TEUCHOS_MAGISTRATE}")
    ADD_SUBDIRECTORY(magistrate)
  ELSE()
    IF(${PACKAGE_NAME}_DOWNLOAD_SOURCE)
      FETCHCONTENT_DECLARE(
        download_magistrate
        GIT_REPOSITORY https://github.com/DARMA-tasking/magistrate
        GIT_TAG db689e717fc0b47e89c554afd87199dbb8486061
      )
      FETCHCONTENT_MAKEAVAILABLE(download_magistrate)
    ELSE()
      MESSAGE(FATAL_ERROR "ERROR: magistrate was enalbed but could not find magistrate source directory. Please check out the repo under the Trilinos/packages/magistrate source directory or set ${PACKAGE_NAME}_DOWNLOAD_SOURCE=ON to automatically download.")
    ENDIF()
  ENDIF()
ENDIF()

ADD_LIBRARY(magistrate::all_libs ALIAS ${MAGISTRATE_LIBRARY})

TRIBITS_PACKAGE_POSTPROCESS()
