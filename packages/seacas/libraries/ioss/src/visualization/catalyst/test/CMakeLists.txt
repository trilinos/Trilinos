CMAKE_MINIMUM_REQUIRED(VERSION 3.8.0)

INCLUDE(CTest)

  set(PhactoriTestDirectory
    ${CMAKE_CURRENT_BINARY_DIR}/PhactoriTests)
  set(PhactoriUnitTestsDirectory
    ${PhactoriTestDirectory}/UnitTests)
  file(MAKE_DIRECTORY ${PhactoriUnitTestsDirectory})

FUNCTION(addCatalystTest TEST_NAME MPI_RANKS_LIST TEST_EXECUTABLE
         EXTRA_TEST_ARGS)
    FOREACH(PROC ${MPI_RANKS_LIST})
        set(localWorkingDir ${CMAKE_CURRENT_BINARY_DIR}/testdirs/${TEST_NAME}_NP_${PROC})
        file(MAKE_DIRECTORY ${localWorkingDir})
        ADD_TEST(NAME ${TEST_NAME}_NP_${PROC}
             COMMAND ${MPIEXEC_EXECUTABLE}
             ${MPIEXEC_NUMPROC_FLAG} ${PROC}
             ${TEST_EXECUTABLE} ${EXTRA_TEST_ARGS}
             WORKING_DIRECTORY ${localWorkingDir}
        )
    ENDFOREACH()
ENDFUNCTION()

ADD_SUBDIRECTORY(fixture)
ADD_SUBDIRECTORY(logging)
ADD_SUBDIRECTORY(parser)
ADD_SUBDIRECTORY(phactori)
ADD_SUBDIRECTORY(pipeline)
