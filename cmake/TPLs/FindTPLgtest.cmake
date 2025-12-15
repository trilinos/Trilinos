set(REQUIRED_HEADERS gtest.h)
set(REQUIRED_LIBS_NAMES "gtest googletest")
set(IMPORTED_TARGETS_FOR_ALL_LIBS GTest::gtest)

tribits_tpl_allow_pre_find_package(gtest GTEST_ALLOW_PREFIND)

set(gtest_ALLOW_DOWNLOAD "ON" CACHE BOOL "Allow CMake to download gtest")

if (GTEST_ALLOW_PREFIND)
  message("-- Using find_package(GTest ...) ...")
  find_package(GTest)
  if (GTest_FOUND)
    message("-- Found GTest_DIR='${GTest_DIR}'")
    message("-- Generating gtest::all_libs and gtestConfig.cmake")
    tribits_extpkg_create_imported_all_libs_target_and_config_file(gtest
      INNER_FIND_PACKAGE_NAME GTest
      IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS})
  else()

    if (gtest_ALLOW_DOWNLOAD)

      include(FetchContent)
      FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        52eb8108c5bdec04579160ae17225d66034bd723 # release-1.17.0
        FIND_PACKAGE_ARGS
      )
      FetchContent_MakeAvailable(googletest)
      find_package(googletest REQUIRED)
      if (googletest_FOUND)
        message("-- Found GTest_DIR='${GTest_DIR}'")
        message("-- Generating gtest::all_libs and gtestConfig.cmake")
        tribits_extpkg_create_imported_all_libs_target_and_config_file(gtest
          INNER_FIND_PACKAGE_NAME googletest
          IMPORTED_TARGETS_FOR_ALL_LIBS  ${IMPORTED_TARGETS_FOR_ALL_LIBS})
      endif()
    endif()
  endif()
endif()

if (NOT TARGET gtest::all_libs)

  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES(gtest
    REQUIRED_HEADERS ${REQUIRED_HEADERS}
    REQUIRED_LIBS_NAMES ${REQUIRED_LIBS_NAMES}
  )

endif()
