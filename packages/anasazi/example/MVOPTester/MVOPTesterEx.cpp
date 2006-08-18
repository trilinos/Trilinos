#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziBasicOutputManager.hpp"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  // If we aren't using MPI, then setup a serial communicator.
  Epetra_SerialComm Comm;
#endif

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::BasicOutputManager<double>() );
  MyOM->setVerbosity( Anasazi::Warnings );
  MyOM->stream(Anasazi::Errors) << Anasazi::Anasazi_Version() << endl << endl;

  // Number of global elements
  int dim = 100;
  int blockSize = 5;

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Epetra_Map Map(dim, 0, Comm);
  
  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map.NumMyElements();
  int * MyGlobalElements = new int[NumMyElements];
  Map.MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  int * NumNz = new int[NumMyElements];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)
  for (int i=0; i<NumMyElements; i++) {
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == dim-1) {
      NumNz[i] = 2;
    }
    else {
      NumNz[i] = 3;
    }
  }

  // Create an Epetra_Matrix
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A 
    = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, NumNz) );
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1
  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
  double two = 2.0;
  int NumEntries;
  for (int i=0; i<NumMyElements; i++) {
    if (MyGlobalElements[i]==0) {
      Indices[0] = 1;
      NumEntries = 1;
    }
    else if (MyGlobalElements[i] == dim-1) {
      Indices[0] = dim-2;
      NumEntries = 1;
    }
    else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
    }
    int ierr = A->InsertGlobalValues(MyGlobalElements[i],NumEntries,Values,Indices);
    assert(ierr==0);
    // Put in the diagonal entry
    ierr = A->InsertGlobalValues(MyGlobalElements[i],1,&two,MyGlobalElements+i);
    assert(ierr==0);
  }
   
  // Finish building the epetra matrix A
  int ierr = A->FillComplete();
  assert(ierr==0);

  // Issue several useful typedefs;
  // The MultiVecTraits class is for defining ....
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note that this needs to have the same number of columns as the blocksize.
  Teuchos::RefCountPtr<MV> ivec 
    = Teuchos::rcp( new MV(Map, blockSize) );
  ivec->Random();

  // test the multivector class
  bool testret;
  testret = Anasazi::TestMultiVecTraits<double,MV>(MyOM,ivec);
  if (testret) {
    MyOM->print(Anasazi::Warnings,"*** PASSED TestMultiVecTraits() ***\n");
  }
  else {
    MyOM->print(Anasazi::Warnings,"*** FAILED TestMultiVecTraits() ***\n");
  }


  // test the operator class
  testret = Anasazi::TestOperatorTraits<double,MV,OP>(MyOM,ivec,A);
  if (testret) {
    MyOM->print(Anasazi::Warnings,"*** PASSED TestOperatorTraits() ***\n");
  }
  else {
    MyOM->print(Anasazi::Warnings,"*** FAILED TestOperatorTraits() ***\n");
  }

  // Release all objects
  delete [] NumNz;
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
