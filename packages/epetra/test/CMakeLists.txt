INCLUDE(${CMAKE_MODULE_PATH}/Epetra_Add_Test.cmake)

# Syntax for adding a test:
#   EPETRA_ADD_TEST(<TEST_DIR_NAME> SOURCES <SOURCE1> <SOURCE2> ... ARGS <ARG1> <ARG2> ...)

# The default arguments include "-v", as specified in the macro EPETRA_ADD_TEST

EPETRA_ADD_TEST(BlockMap SOURCES cxx_main.cpp checkmap.cpp ARGS -v)

EPETRA_ADD_TEST(BasicPerfTest SOURCES cxx_main.cpp ARGS 16 12 1 1 25 -v INSTALL)

EPETRA_ADD_TEST(Comm SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(CrsGraph SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(CrsMatrix SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(RowMatrix SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(CrsRectMatrix SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(Directory SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(FECrsGraph SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)

EPETRA_ADD_TEST(FECrsMatrix SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)

EPETRA_ADD_TEST(FEVbrMatrix SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)

EPETRA_ADD_TEST(FEVector SOURCES cxx_main.cpp ExecuteTestProblems.cpp ARGS -v)

EPETRA_ADD_TEST(ImportExport SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(IntSerialDense SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(Map SOURCES cxx_main.cpp checkmap.cpp ARGS -v)

# This test was commented out in the AUTOTOOLS configuration:
EPETRA_ADD_TEST(MapColoring SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(MultiVector SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)

EPETRA_ADD_TEST(Object SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(RowMatrixTransposer SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(SerialDense SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(SerialSpdDense SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(VbrMatrix SOURCES cxx_main.cpp ARGS -v)

EPETRA_ADD_TEST(Vector SOURCES cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp ARGS -v)
