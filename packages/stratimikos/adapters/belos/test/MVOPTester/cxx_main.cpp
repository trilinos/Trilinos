// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
//  This test uses the MVOPTester.hpp functions to test the Belos adapters
//  to Epetra and Thyra.
//

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosEpetraAdapter.hpp"

#ifdef HAVE_EPETRA_THYRA
#include "BelosThyraAdapter.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif



int main(int argc, char *argv[])
{

  using Teuchos::rcp_implicit_cast;

  int i, ierr, gerr;
  gerr = 0;
#ifdef NDEBUG
  (void)ierr;  // Eliminate unused warning
#endif

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Teuchos::RCP<Epetra_MpiComm> Comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  // If we aren't using MPI, then setup a serial communicator.
  Teuchos::RCP<Epetra_SerialComm> Comm = Teuchos::rcp( new Epetra_SerialComm() );
#endif


   // number of global elements
  int dim = 100;
  int blockSize = 3;

  // PID info
  int MyPID = Comm->MyPID();
  bool verbose = 0;

  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Teuchos::RCP<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );

  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map->NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map->MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer std::vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor
  std::vector<int> NumNz(NumMyElements);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)
  for (i=0; i<NumMyElements; i++) {
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == dim-1) {
      NumNz[i] = 2;
    }
    else {
      NumNz[i] = 3;
    }
  }

  // Create an Epetra_Matrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *Map, &NumNz[0]) );

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1
  std::vector<double> Values(2);
  Values[0] = -1.0; Values[1] = -1.0;
  std::vector<int> Indices(2);
  double two = 2.0;
  int NumEntries;
  for (i=0; i<NumMyElements; i++) {
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
    ierr = A->InsertGlobalValues(MyGlobalElements[i],NumEntries,&Values[0],&Indices[0]);
    assert(ierr==0);
    // Put in the diagonal entry
    ierr = A->InsertGlobalValues(MyGlobalElements[i],1,&two,&MyGlobalElements[i]);
    assert(ierr==0);
  }

  // Finish building the epetra matrix A
  ierr = A->FillComplete();
  assert(ierr==0);

  // Create an Belos::EpetraOp from this Epetra_CrsMatrix
  Teuchos::RCP<Belos::EpetraOp> op = Teuchos::rcp(new Belos::EpetraOp(A));

  // Issue several useful typedefs;
  // typedef Belos::MultiVec<double> EMV; // unused
  // typedef Belos::Operator<double> EOP; // unused

  // Create an Epetra_MultiVector for an initial std::vector to start the solver.
  // Note that this needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Belos::EpetraMultiVec> ivec = Teuchos::rcp( new Belos::EpetraMultiVec(*Map, blockSize) );
  ivec->Random();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Belos::OutputManager<double> > MyOM = Teuchos::rcp( new Belos::OutputManager<double>( MyPID ) );
  if (verbose) {
    MyOM->setVerbosity( Belos::Errors + Belos::Warnings );
  }

#ifdef HAVE_EPETRA_THYRA
  typedef Thyra::MultiVectorBase<double> TMVB;
  typedef Thyra::LinearOpBase<double>    TLOB;
  // create thyra objects from the epetra objects

  // first, a Thyra::VectorSpaceBase
  Teuchos::RCP<const Thyra::VectorSpaceBase<double> > epetra_vs =
    Thyra::create_VectorSpace(Map);

  // then, a MultiVectorBase (from the Epetra_MultiVector)
  Teuchos::RCP<Thyra::MultiVectorBase<double> > thyra_ivec =
    Thyra::create_MultiVector(rcp_implicit_cast<Epetra_MultiVector>(ivec),epetra_vs);

  // then, a LinearOpBase (from the Epetra_CrsMatrix)
  Teuchos::RCP<Thyra::LinearOpBase<double> > thyra_op =
    Teuchos::rcp( new Thyra::EpetraLinearOp(A) );


  // test the Thyra adapter multivector
  ierr = Belos::TestMultiVecTraits<double,TMVB>(MyOM,thyra_ivec);
  gerr |= ierr;
  switch (ierr) {
  case Belos::Ok:
    if ( verbose && MyPID==0 ) {
      std::cout << "*** ThyraAdapter PASSED TestMultiVecTraits()" << std::endl;
    }
    break;
  case Belos::Error:
    if ( verbose && MyPID==0 ) {
      std::cout << "*** ThyraAdapter FAILED TestMultiVecTraits() ***"
           << std::endl << std::endl;
    }
    break;
  }

  // test the Thyra adapter operator
  ierr = Belos::TestOperatorTraits<double,TMVB,TLOB>(MyOM,thyra_ivec,thyra_op);
  gerr |= ierr;
  switch (ierr) {
  case Belos::Ok:
    if ( verbose && MyPID==0 ) {
      std::cout << "*** ThyraAdapter PASSED TestOperatorTraits()" << std::endl;
    }
    break;
  case Belos::Error:
    if ( verbose && MyPID==0 ) {
      std::cout << "*** ThyraAdapter FAILED TestOperatorTraits() ***"
           << std::endl << std::endl;
    }
    break;
  }
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (gerr) {
    if (verbose && MyPID==0)
      std::cout << "End Result: TEST FAILED" << std::endl;
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;

}
