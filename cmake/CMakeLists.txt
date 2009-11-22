
# Label all of these tests as TrilinosFramework tests!
SET(PACKAGE_NAME TrilinosFramework)


IF (Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS)
  ADD_SUBDIRECTORY(DependencyUnitTests)
ENDIF()

IF (Trilinos_ENABLE_TESTING_UNIT_TESTS)
  ADD_SUBDIRECTORY(TestingUnitTests)
ENDIF()

IF (Trilinos_ENABLE_PYTHON_UNIT_TESTS)
  ADD_SUBDIRECTORY(python)
ENDIF()
