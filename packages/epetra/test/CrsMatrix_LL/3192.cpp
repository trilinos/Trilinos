#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {

#ifdef EPETRA_MPI
void
gathervPrintMpi (std::ostream& out, const std::string& str, MPI_Comm comm)
{
  const int rootRank = 0;
  int myRank = 0;
  MPI_Comm_rank (comm, &myRank);

  if (myRank == rootRank) {
    out << str;
  }
  int numProcs;
  MPI_Comm_size (comm, &numProcs);

  const int sizeTag = 156;
  const int msgTag = 418;
  MPI_Status status;
  for (int sendRank = 1; sendRank < numProcs; ++sendRank) {
    if (myRank == sendRank) {
      int msgSize = static_cast<int> (str.size ());
      MPI_Send (&msgSize, 1, MPI_INT, rootRank, sizeTag, comm);
      if (msgSize != 0) {
        MPI_Send (str.c_str (), msgSize+1, MPI_CHAR, rootRank, msgTag, comm);
      }
    }
    else if (myRank == rootRank) {
      int msgSize = 0;
      MPI_Recv (&msgSize, 1, MPI_INT, sendRank, sizeTag, comm, &status);
      if (msgSize != 0) { // msgSize doesn't count '\0'
        std::vector<char> msgBuf (msgSize + std::size_t (1));
        MPI_Recv (msgBuf.data (), msgSize+1, MPI_CHAR,
                  sendRank, msgTag, comm, &status);
        std::string msg (msgBuf.data ());
        out << msg;
      }
    }
  }
}
#endif // EPETRA_MPI

void
gathervPrint (std::ostream& out,
              const std::string& str,
              const Epetra_Comm&
#ifdef EPETRA_MPI
              comm
#endif // EPETRA_MPI
              )
{
#ifdef EPETRA_MPI
  const Epetra_MpiComm* mpiComm = dynamic_cast<const Epetra_MpiComm*> (&comm);
  if (mpiComm != nullptr) {
    gathervPrintMpi (out, str, mpiComm->Comm ());
    return;
  }
#endif // EPETRA_MPI
  out << str;
}

}

template<class GlobalOrdinalType>
int
testEpetra (std::ostream& out,
            const char globalOrdinalTypeName[],
            const Epetra_Comm& comm)
{
  using std::endl;
  using GO = GlobalOrdinalType;
  out << "Test Epetra with global ordinal type " << globalOrdinalTypeName
      << " (sizeof is " << sizeof (GlobalOrdinalType) << ")" << endl;
  const GO localRows = 5000;
  const int myRank = comm.MyPID ();
  int lclSuccess = 1;
  int gblSuccess = 1;
  int epetraErrCode = 0;

  if (myRank == 0) {
    out << "Create row Map" << endl;
  }
  std::unique_ptr<Epetra_Map> rowMap;
  std::ostringstream errMsg;
  try {
    rowMap = std::unique_ptr<Epetra_Map> (new Epetra_Map (localRows, 0, comm));
  }
  catch (std::exception& e) {
    errMsg << "Proc " << myRank << ": Epetra_Map constructor threw an "
      "exception: " << e.what () << endl;
    lclSuccess = 0;
  }
  catch (const char* errStr) {
    errMsg << "Proc " << myRank << ": Epetra_Map constructor threw "
      "const char*: " << errStr << endl;
    lclSuccess = 0;
  }
  catch (...) {
    errMsg << "Proc " << myRank << ": Epetra_Map constructor threw an "
      "exception not a subclass of std::exception" << endl;
    lclSuccess = 0;
  }
  comm.MinAll (&lclSuccess, &gblSuccess, 1);
  if (gblSuccess != 1) {
    out << "Epetra_Map constructor threw an exception!" << endl;
    gathervPrint (out, errMsg.str (), comm);
    return gblSuccess;
  }
  const int gblNumRows = rowMap->NumGlobalElements ();
  const int lclNumRows = rowMap->NumMyElements ();
  if (myRank == 0) {
    out << "Global number of rows: " << gblNumRows << endl
        << "Local number of rows: " << lclNumRows << endl;
  }

  if (myRank == 0) {
    out << "Create CrsMatrix" << endl;
  }
  std::unique_ptr<Epetra_CrsMatrix> A;
  try {
    A = std::unique_ptr<Epetra_CrsMatrix> (new Epetra_CrsMatrix (Epetra_DataAccess::Copy, *rowMap, 3));
  }
  catch (std::exception& e) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix constructor threw an "
      "exception: " << e.what () << endl;
    lclSuccess = 0;
  }
  catch (const char* errStr) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix constructor threw "
      "const char*: " << errStr << endl;
    lclSuccess = 0;
  }
  catch (...) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix constructor threw an "
      "exception not a subclass of std::exception" << endl;
    lclSuccess = 0;
  }
  comm.MinAll (&lclSuccess, &gblSuccess, 1);
  if (gblSuccess != 1) {
    out << "Epetra_CrsMatrix constructor failed!" << endl;
    gathervPrint (out, errMsg.str (), comm);
    return gblSuccess;
  }

  if (myRank == 0) {
    out << "Insert values into CrsMatrix" << endl;
  }
  double tile[3] = {-2.0, 4.0, -2.0};
  // Per #3192, we deliberately use int here, regardless of GlobalOrdinalType.
  int firstCols[] = {0, 1};
  try {
    epetraErrCode = A->InsertGlobalValues (0, 2, tile + 1, firstCols);
    if (epetraErrCode != 0) {
      errMsg << "Proc " << myRank << ": First InsertGlobalValues failed "
        "with error code " << epetraErrCode << endl;
      lclSuccess = 0;
    }

    epetraErrCode = 0;
    for (int row = 1; row < lclNumRows - 1; ++row) {
      int cols[] = {row - 1, row, row + 1};
      const int err = A->InsertGlobalValues (row, 3, tile, cols);
      if (err != 0) {
        lclSuccess = 0;
        if (epetraErrCode == 0) {
          epetraErrCode = err;
        }
      }
    }

    if (epetraErrCode != 0) {
      errMsg << "Proc " << myRank << ": At least one middle InsertGlobalValues "
        "call failed with nonzero error code " << epetraErrCode << endl;
    }
    int lastCols[] = {lclNumRows - 2, lclNumRows - 1};
    epetraErrCode = A->InsertGlobalValues (lclNumRows - 1, 2, tile, lastCols);
    if (epetraErrCode != 0) {
      errMsg << "Proc " << myRank << ": Last InsertGlobalValues failed "
        "with error code " << epetraErrCode << endl;
      lclSuccess = 0;
    }
  }
  catch (std::exception& e) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::InsertGlobalValues "
      "threw an exception: " << e.what () << endl;
    lclSuccess = 0;
  }
  catch (const char* errStr) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::InsertGlobalValues "
      "threw const char*: " << errStr << endl;
    lclSuccess = 0;
  }
  catch (...) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::InsertGlobalValues "
      "threw an exception not a subclass of std::exception" << endl;
    lclSuccess = 0;
  }
  comm.MinAll (&lclSuccess, &gblSuccess, 1);
  if (gblSuccess != 1) {
    out << "Epetra_CrsMatrix::InsertGlobalValues failed!" << endl;
    gathervPrint (out, errMsg.str (), comm);
    return gblSuccess;
  }

  if (myRank == 0) {
    out << "Call FillComplete on the matrix" << endl;
  }
  try {
    epetraErrCode = A->FillComplete ();
    if (epetraErrCode != 0) {
      errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::FillComplete "
        "returned a nonzero error code " << epetraErrCode << endl;
      lclSuccess = 0;
    }
  }
  catch (std::exception& e) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::FillComplete "
      "threw an exception: " << e.what () << endl;
    lclSuccess = 0;
  }
  catch (const char* errStr) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::FillComplete threw "
      "const char*: " << errStr << endl;
    lclSuccess = 0;
  }
  catch (...) {
    errMsg << "Proc " << myRank << ": Epetra_CrsMatrix::FillComplete "
      "threw an exception not a subclass of std::exception" << endl;
    lclSuccess = 0;
  }
  comm.MinAll (&lclSuccess, &gblSuccess, 1);
  if (gblSuccess != 1) {
    out << "Epetra_CrsMatrix::FillComplete failed!" << endl;
    gathervPrint (out, errMsg.str (), comm);
    return gblSuccess;
  }

  return gblSuccess;
}

int
main (int argc, char *argv[])
{
#ifdef EPETRA_MPI
  MPI_Init (&argc, &argv);
  MPI_Comm mpiComm = MPI_COMM_WORLD;
  Epetra_MpiComm comm (mpiComm);
#else
  Epetra_SerialComm comm;
#endif // EPETRA_MPI

  const int gbl32bitResult =
    testEpetra<int> (std::cerr, "int", comm);
  const int gbl64bitResult =
    testEpetra<long long> (std::cerr, "long long", comm);
  const bool success = gbl32bitResult == 1 && gbl64bitResult == 1;
  if (comm.MyPID () == 0) {
    if (success) {
      std::cout << "Test PASSED" << std::endl;
    }
    else {
      std::cout << "Test FAILED" << std::endl;
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize ();
#endif // EPETRA_MPI

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
