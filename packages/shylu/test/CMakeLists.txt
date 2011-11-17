
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    aztec_driver
    SOURCES shylu_driver.cpp
    DEPLIBS isorropia ifpack aztecoo zoltan belos ml amesos epetra teuchos
    COMM mpi
    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    belos_driver
    SOURCES shylu_belos_driver.cpp
    DEPLIBS isorropia ifpack aztecoo zoltan belos ml amesos epetra teuchos
    COMM mpi
    )

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUDriverFiles
    SOURCE_FILES ShyLU.xml wathen96.mtx wathenSmall.mtx wathenMedium.mtx wathenLarge.mtx
    EXEDEPS aztec_driver belos_driver
    )
