INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    MeshTyingBlocked_SimpleSmoother
#    SOURCES MeshTyingBlocked_SimpleSmoother.cpp
#    COMM serial mpi
#    )

TRIBITS_ADD_EXECUTABLE(
  MyFirstExecutable
  SOURCES main.cpp
  COMM serial mpi
)

TRIBITS_ADD_TEST(
  MyFirstExecutable
  NAME "TestMyFirstExecutable"
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 4
)

#  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tutorial_ConcatenateFragmentsCp
#    SOURCE_FILES
#      concatenate_fragment_files.py
#    )
#
#  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tutorial_SampleTestCp
#    SOURCE_FILES
#      main_sample_test.cpp.fragment.01
#      main_sample_test.cpp.fragment.02
#    )
#
#  TRIBITS_ADD_ADVANCED_TEST(TutorialSampleTest
#    TEST_0
#      CMND ${PYTHON_EXECUTABLE}
#      ARGS concatenate_fragment_files.py main_sample_test.cpp 2 main_sample_test.cpp.fragment.01 main_sample_test.cpp.fragment.02
#    )
