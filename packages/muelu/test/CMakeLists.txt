ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Tpetra
  ${PACKAGE_NAME}_ENABLE_Ifpack2
  ${PACKAGE_NAME}_ENABLE_Amesos2
  Xpetra_INT_LONG_LONG
  )

ADD_SUBDIRECTORY(unit_tests)

IF (${PACKAGE_NAME}_ENABLE_Kokkos_Refactor)
    ADD_SUBDIRECTORY(unit_tests_kokkos)
    ADD_SUBDIRECTORY(perf_tests_kokkos)
ENDIF()

ADD_SUBDIRECTORY(scaling)

if (NOT MUELU_VALGRIND_TESTING)
    ADD_SUBDIRECTORY(interface)
ELSE()
    MESSAGE(STATUS "Skipping ParameterListInterpreter.")
ENDIF()

ADD_SUBDIRECTORY(blockedtransfer)

IF(Xpetra_INT_LONG_LONG AND HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
    ADD_SUBDIRECTORY(longlong)
ENDIF()

ADD_SUBDIRECTORY(factories)
ADD_SUBDIRECTORY(maxwell)
ADD_SUBDIRECTORY(navierstokes)

if (NOT MUELU_VALGRIND_TESTING)
    # JJH 4Feb2014:  This is failing the nightly MPI valgrind test.  It runs fine
    # interactively.  If/when I figure out the cause, I'll reenable this during
    # memory testing.
    ADD_SUBDIRECTORY(paramlist)
ENDIF()

ADD_SUBDIRECTORY(permutation)
ADD_SUBDIRECTORY(profiling)
ADD_SUBDIRECTORY(simple1D)

IF(NOT ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor)
    # AggregationExportFactory can't work with
    # Aggregates_kokkos (yet)
    ADD_SUBDIRECTORY(viz)
ENDIF()

ADD_SUBDIRECTORY(toggletransfer)
ADD_SUBDIRECTORY(perf_tests)

if (NOT MUELU_VALGRIND_TESTING)
    ADD_SUBDIRECTORY(vardofpernode)
ELSE()
    MESSAGE(STATUS "Skipping vardofpernode test due to valgrind.")
ENDIF()

if (NOT MUELU_VALGRIND_TESTING)
    ADD_SUBDIRECTORY(convergence)
    IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
        ADD_SUBDIRECTORY(helmholtz)
    ENDIF()
ELSE()
    MESSAGE(STATUS "Skipping Helmholtz due to valgrind testing.")
ENDIF()

#ADD_SUBDIRECTORY(mhdvanka)

ADD_SUBDIRECTORY(structured)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)

SET(SOURCES "")
