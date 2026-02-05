TRIBITS_PACKAGE(magistrate)

include(FetchContent)

if(${PROJECT_NAME}_ENABLE_magistrate)

  # Disable target exports in magistrate so that tribits can export the
  # targets as a trilinos package
  set(magistrate_enable_export_of_targets OFF)

  if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/magistrate")
    set(MAGISTRATE_REPO_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/magistrate")
    message(STATUS "magistrate package found in Trilinos source")
    message(STATUS "magistrate source dir = ${MAGISTRATE_REPO_SOURCE_DIR}")
    add_subdirectory(magistrate)
  else()
    message(STATUS "Downloading the magistrate package with FetchContent.")
    fetchcontent_declare(magistrate_repo
      GIT_REPOSITORY https://github.com/darma-tasking/magistrate.git
      GIT_TAG develop
    )
    fetchcontent_makeavailable(magistrate_repo)
    fetchcontent_getproperties(magistrate_repo SOURCE_DIR MAGISTRATE_REPO_SOURCE_DIR)
    message(STATUS "magistrate source dir = ${MAGISTRATE_REPO_SOURCE_DIR}")
  endif()

  add_library(magistrate::all_libs ALIAS ${MAGISTRATE_LIBRARY})

endif()

TRIBITS_PACKAGE_POSTPROCESS()
