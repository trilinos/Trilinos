// includes for the test
#include <iostream>
#include <Epetra_Version.h>
#include <Epetra_SerialComm.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <Epetra_Time.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <stdlib.h>
#include <vector>

#include "TestUtils.h"
#include "TPICrsMatrix.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif

using std::endl;
using std::cout;
using std::string;
using std::setw;
using std::setprecision;
using std::fixed;
using std::scientific;
using std::vector;

#define PRINT() \
{ \
    if (rank == 0) cout << "Usage: " << argv[0] << " [matrix filename] [thread args] [numTrials]" << endl \
                        << "where [thread args] takes one of the following forms:" << endl \
                        << " b:e+i    e.g., 1:5+1 runs 1,2,3,4,5 threads" << endl \
                        << " b:e*m    e.g., 1:2:16 runs 1,2,4,8,16 threads" << endl \
                        << " n[,n]*   e.g., 1,3,9,12 runs 1,3,9,12 threads" << endl; \
}

#ifdef EPETRA_MPI
  #define PRINT_AND_EXIT() \
  { \
     PRINT(); \
     MPI_Finalize(); \
     return -1; \
  }
#else 
  #define PRINT_AND_EXIT() \
  { \
     PRINT(); \
     return -1; \
  }
#endif

void getStats(Epetra_CrsMatrix *&L, const char *fn, int *&NNZperRow) {
  cout << endl << "*** Matrix statistics" << endl;
  // FINISH: gather/output/return statitics
  Epetra_SerialComm Comm;
  if (L == NULL) {
    // read from file
    int ierr = EpetraExt::MatrixMarketFileToCrsMatrix(fn,Comm,L,false,false); // transpose=false, verbose=false
    if (ierr) {
      cout << "EpetraExt::MatrixMarketFIleToCrsMatrix returned error " << ierr << "." << endl;
      throw std::runtime_error("file error");
    }
  }
  const int NNZ = L->NumMyNonzeros();
  const int NLRs = L->NumMyRows();
  NNZperRow = new int[NLRs];
  double stddev = 0;
  double mean = (double)(NNZ) / (double)(NLRs);
  for (int i=0; i<NLRs; ++i) {
    NNZperRow[i] = L->NumMyEntries(i);
    double tmp = (NNZperRow[i] - mean);
    stddev += (tmp*tmp);
  }
  stddev = sqrt(stddev / NLRs);
  cout << "Number of rows: " << NLRs << endl;
  cout << "Number of non-zeros: " << NNZ << endl;
  cout << "Mean number of non-zeros per row: " << fixed << setprecision(1) << mean << endl;
  cout << "Std dev number of non-zeros per row: " << fixed << setprecision(2) << stddev << endl;
}


void runEpetraTest(int rank, int size, Epetra_CrsMatrix *L, int numTrials) {
  if (rank == 0) cout << endl << "*** Epetra_CrsMatrix over " << numTrials << " trials" << endl;
  Epetra_Time timer(L->Comm());
  // test LevelSolver
  Epetra_Vector x(L->RowMap(),false); 
  Epetra_Vector Lx(x);
  x.Random();
  Lx.Random();
  double time = 0;
  timer.ResetStartTime();
  for (int t=0; t<numTrials; ++t) {
    L->Apply(Lx,x);
  }
  time = timer.ElapsedTime();
  if (rank == 0) {
    cout << setw(20) << "Epetra::Apply(), " << setw(2) << size << " processes, " 
         << setprecision(2) << scientific << setw(11) << time            << " total, " 
         << setprecision(2) << scientific << setw(11) << time /numTrials << " average" 
         << endl;
  }
}

// TODO: separate this; no need to rebuild the matrix for every test run
void runTPITest(int numThreads, Epetra_CrsMatrix *&L, const int *NNZperRow, int numTrials) {
  Epetra_SerialComm Comm;
  Epetra_Time timer(Comm);
  // get NNZ per row info needed for TPICrsMatrix allocation
  const int NLR = L->RowMap().NumMyPoints();
  TPICrsMatrix Lpar(L->RowMap(),L->ColMap(),NNZperRow);
  // put matrix data into TPICrsMatrix
  for (int r=0; r<NLR; ++r) {
    double *vals;
    int    *inds;
    int     numentries;
    L->ExtractMyRowView(r,numentries,vals,inds);
    assert( Lpar.InsertMyValues(r,numentries,vals,inds) == 0);
  }
  assert( Lpar.FillComplete(L->OperatorDomainMap(),L->OperatorRangeMap()) == 0);
  // test LevelSolver
  Epetra_Vector x(L->RowMap(),false); 
  Epetra_Vector Lx(x);
  x.Random();
  Lx.Random();
  int ret = TPI_Init(numThreads);
  if (ret != numThreads) { cout << "error " << ret << " from TPI_Init" << endl; throw std::runtime_error("TPI error"); }
  double time = 0;
  timer.ResetStartTime();
  for (int t=0; t<numTrials; ++t) {
    Lpar.Apply(x,Lx);
  }
  time = timer.ElapsedTime();
  cout << endl << setw(20) << "TPICrs::Apply(), " << setw(2) << numThreads << " threads,   " 
    << setprecision(2) << scientific << setw(11) << time            << " total, " 
    << setprecision(2) << scientific << setw(11) << time /numTrials << " average" 
    << endl;
  ret = TPI_Finalize();
  if (ret) { cout << "error " << ret << " from TPI_Finalize" << endl; throw std::runtime_error("TPI error"); }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main() routines for the test
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  int rank, size;
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  rank = 0;
  size = 1;
  Epetra_SerialComm Comm;
#endif
  Epetra_Time time(Comm);
  // if (rank == 0) cout << Epetra_Version() << endl << endl;

  int numTrials = 100;
  string fn;
  vector<int> threads;
  if (argc > 1) {
    fn = argv[1];
    if (argc > 2) {
      if ( parseNumThreads(argv[2],threads) != 0 ) {
        PRINT_AND_EXIT();
      }
    }
    if (argc > 3) {
      if (sscanf(argv[3],"%d",&numTrials) != 1) PRINT_AND_EXIT();
    }
  }
  else {
    PRINT_AND_EXIT();
  }

  // read matrix data into Epetra_CrsMatrix object
  // not only does this let us use the EpetraExt functionality, but we can use the 
  // resulting matrix to help us test our TPICrsMatrix
  Epetra_CrsMatrix *L;
  int ierr =
    EpetraExt::MatrixMarketFileToCrsMatrix(fn.c_str(),Comm,L,false,false); // transpose=false, verbose=false
  if (ierr) {
    if (rank == 0) cout << "EpetraExt::MatrixMarketFIleToCrsMatrix returned error " << ierr << "." << endl;
    return -1;
  }

  runEpetraTest(rank,size,L,numTrials);

  if (rank == 0 && !threads.empty()) {
    int *NNZperRow;
    if (size > 1) {
      // getStats will create a new Epetra_CrsMatrix, which wlll be deleted below
      delete L;
      L = NULL;
    }
    getStats(L,fn.c_str(),NNZperRow);
    cout << "*** TPICrsMatrix over " << numTrials << " trials" << endl;
    for (vector<int>::iterator nt=threads.begin(); nt != threads.end(); ++nt) {
      runTPITest(*nt,L,NNZperRow,numTrials);
    }
    delete NNZperRow;
  }
  
  // delete manually allocated matrices
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  delete L;
  return 0;
}
