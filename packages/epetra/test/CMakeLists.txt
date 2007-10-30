INCLUDE(${CMAKE_MODULE_PATH}/Epetra_Add_Test.cmake)

# Syntax for adding a test:
#   To specify the arguments to run the test with:
#     SET(EXTRA_TEST_ARGS <P1> <P2> ...)
#   To add the test:
#     EPETRA_ADD_TEST(<TEST_DIR_NAME> <TEST_NAME> <SOURCE1> <SOURCE2> ...)

# The default arguments include "-v", as specified in the macro EPETRA_ADD_TEST

EPETRA_ADD_TEST(BlockMap BlockMap cxx_main.cpp checkmap.cpp)

SET(TEST_ARGS 16 12 1 1 25)
EPETRA_ADD_TEST(BasicPerfTest BasicPerfTest cxx_main.cpp)

EPETRA_ADD_TEST(Comm Comm cxx_main.cpp)

EPETRA_ADD_TEST(CrsGraph CrsGraph cxx_main.cpp)

EPETRA_ADD_TEST(CrsMatrix CrsMatrix cxx_main.cpp)

EPETRA_ADD_TEST(RowMatrix RowMatrix cxx_main.cpp)

EPETRA_ADD_TEST(CrsRectMatrix CrsRectMatrix cxx_main.cpp)

EPETRA_ADD_TEST(Directory Directory cxx_main.cpp)

EPETRA_ADD_TEST(FECrsGraph FECrsGraph cxx_main.cpp ExecuteTestProblems.cpp)

EPETRA_ADD_TEST(FECrsMatrix FECrsMatrix cxx_main.cpp ExecuteTestProblems.cpp)

EPETRA_ADD_TEST(FEVbrMatrix FEVbrMatrix cxx_main.cpp ExecuteTestProblems.cpp)

EPETRA_ADD_TEST(FEVector FEVector cxx_main.cpp ExecuteTestProblems.cpp)

EPETRA_ADD_TEST(ImportExport ImportExport cxx_main.cpp)

EPETRA_ADD_TEST(IntSerialDense IntSerialDense cxx_main.cpp)

EPETRA_ADD_TEST(Map Map cxx_main.cpp checkmap.cpp)

# This test was commented out in the AUTOTOOLS configuration:
EPETRA_ADD_TEST(MapColoring MapColoring cxx_main.cpp)

EPETRA_ADD_TEST(MultiVector MultiVector cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp)

EPETRA_ADD_TEST(Object Object cxx_main.cpp)

EPETRA_ADD_TEST(RowMatrixTransposer RowMatrixTransposer cxx_main.cpp)

EPETRA_ADD_TEST(SerialDense SerialDense cxx_main.cpp)

EPETRA_ADD_TEST(SerialSpdDense SerialSpdDense cxx_main.cpp)

EPETRA_ADD_TEST(VbrMatrix VbrMatrix cxx_main.cpp)

EPETRA_ADD_TEST(Vector Vector cxx_main.cpp BuildTestProblems.cpp ExecuteTestProblems.cpp)
