
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    aztec_driver
#    SOURCES shylu_driver.cpp
#    COMM mpi
#    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    belos_driver
    SOURCES shylu_belos_driver.cpp
    COMM mpi
    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    iqr_driver
    SOURCES shylu_iqr_driver.cpp
    COMM mpi
    NUM_MPI_PROCS 4
    )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUDriverFiles
    SOURCE_FILES ShyLU.xml shylu_iqr_parameters.xml
    shylu_no_iqr_parameters.xml wathen96.mtx wathenSmall.mtx
    shylu_iqr_parameters.seq.xml
#    EXEDEPS aztec_driver belos_driver
    EXEDEPS belos_driver iqr_driver
    )
