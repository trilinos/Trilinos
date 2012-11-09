
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

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUDriverFiles
    SOURCE_FILES ShyLU.xml wathen96.mtx wathenSmall.mtx
#    EXEDEPS aztec_driver belos_driver
    EXEDEPS belos_driver
    )
