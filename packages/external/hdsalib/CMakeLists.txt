TRIBITS_PACKAGE(HdsaLib)

include(FetchContent)

if(${PROJECT_NAME}_ENABLE_HdsaLib)

  # Disable target exports in HdsaLib so that tribits can export the
  # targets as a trilinos package
  set(HdsaLib_enable_export_of_targets OFF)

  if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/HdsaLib")
    set(HDSALIB_REPO_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/HdsaLib")
    message(STATUS "HdsaLib package found in Trilinos source")
    message(STATUS "HdsaLib source dir = ${HDSALIB_REPO_SOURCE_DIR}")
  else()
    message(STATUS "Downloading the HdsaLib package with FetchContent.")
    fetchcontent_declare(HdsaLib_repo
      GIT_REPOSITORY https://github.com/sandialabs/HdsaLib
      GIT_TAG master
    )
    fetchcontent_populate(HdsaLib_repo)
    fetchcontent_getproperties(HdsaLib_repo SOURCE_DIR HDSALIB_REPO_SOURCE_DIR)
    message(STATUS "HdsaLib source dir = ${HDSALIB_REPO_SOURCE_DIR}")
  endif()

  set(HDSALIB_REPO_SOURCE_DIR "${HDSALIB_REPO_SOURCE_DIR}" CACHE PATH "HdsaLib source directory")

endif()

TRIBITS_PACKAGE_POSTPROCESS()
