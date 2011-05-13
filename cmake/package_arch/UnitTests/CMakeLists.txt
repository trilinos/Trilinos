
INCLUDE(PackageAddAdvancedTest)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_MPI)

PACKAGE_ADD_ADVANCED_TEST( TestingFunctionMacro_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      -DPACKAGE_ARCH_SOURCE_DIR="${PACKAGE_SOURCE_DIR}/package_arch"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/TestingFunctionMacro_UnitTests.cmake"
    PASS_REGULAR_EXPRESSION_ALL
      "Final UnitTests Result: num_run = 104"
      "Final UnitTests Result: PASSED"
  )

PACKAGE_ADD_ADVANCED_TEST( CompilerOptions_UnitTests
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TEST_0 CMND ${CMAKE_COMMAND}
    ARGS
      -DPACKAGE_ARCH_SOURCE_DIR="${PACKAGE_SOURCE_DIR}/package_arch"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/CompilerOptions_UnitTests.cmake"
    PASS_REGULAR_EXPRESSION_ALL
      "Final UnitTests Result: num_run = 129"
      "Final UnitTests Result: PASSED"
  )
