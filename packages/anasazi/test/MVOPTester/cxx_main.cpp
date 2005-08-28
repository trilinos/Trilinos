//@HEADER
// ************************************************************************
// 
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
//
//  This test uses the MVOPTester.hpp functions to test the Anasazi adapters
//  to Epetra and Thyra.
//

// #define TEST_THYRA_ADAPTER

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

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziEpetraAdapter.hpp"

#ifdef TEST_THYRA_ADAPTER
#include "AnasaziThyraAdapter.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif

int main(int argc, char *argv[])
{
  int i, ierr, gerr;
  gerr = 0;

#ifdef HAVE_MPI
  // Initialize MPI and setup an Epetra communicator
  MPI_Init(&argc,&argv);
  Teuchos::RefCountPtr<Epetra_MpiComm> Comm = Teuchos::rcp( new Epetra_MpiComm(MPI_COMM_WORLD) );
#else
  // If we aren't using MPI, then setup a serial communicator.
  Teuchos::RefCountPtr<Epetra_SerialComm> Comm = Teuchos::rcp( new Epetra_SerialComm() );
#endif

   // number of global elements
  int dim = 100;
  int blockSize = 5;

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
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( new Epetra_Map(dim, 0, *Comm) );
  
  // Get update list and number of local equations from newly created Map.
  int NumMyElements = Map->NumMyElements();
  int * MyGlobalElements = new int[NumMyElements];
  Map->MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor
  int * NumNz = new int[NumMyElements];

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
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, *Map, NumNz) );
   
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1
  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
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
    ierr = A->InsertGlobalValues(MyGlobalElements[i],NumEntries,Values,Indices);
    assert(ierr==0);
    // Put in the diagonal entry
    ierr = A->InsertGlobalValues(MyGlobalElements[i],1,&two,MyGlobalElements+i);
    assert(ierr==0);
  }
   
  // Finish building the epetra matrix A
  ierr = A->FillComplete();
  assert(ierr==0);

  // Create an Anasazi::EpetraSymOp from this Epetra_CrsMatrix
  Teuchos::RefCountPtr<Anasazi::EpetraSymOp> op = Teuchos::rcp(new Anasazi::EpetraSymOp(A));

  // Issue several useful typedefs;
  typedef Anasazi::MultiVec<double> EMV;
  typedef Anasazi::Operator<double> EOP;

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note that this needs to have the same number of columns as the blocksize.
  Teuchos::RefCountPtr<Anasazi::EpetraMultiVec> ivec = Teuchos::rcp( new Anasazi::EpetraMultiVec(*Map, blockSize) );
  ivec->Random();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
  if (verbose) {
    MyOM->SetVerbosity( Anasazi::Warning );
  }

  // test the Epetra adapter multivector
  ierr = Anasazi::TestMultiVecTraits<double,EMV>(MyOM,ivec);
  gerr |= ierr;
  switch (ierr) {
  case Anasazi::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** EpetraAdapter PASSED TestMultiVecTraits()" << endl;
    }
    break;
  case Anasazi::Failed:
    if ( verbose && MyPID==0 ) {
      cout << "*** EpetraAdapter FAILED TestMultiVecTraits() ***" 
           << endl << endl;
    }
    break;
  }

  // test the Epetra adapter operator 
  ierr = Anasazi::TestOperatorTraits<double,EMV,EOP>(MyOM,ivec,op);
  gerr |= ierr;
  switch (ierr) {
  case Anasazi::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** EpetraAdapter PASSED TestOperatorTraits()" << endl;
    }
    break;
  case Anasazi::Failed:
    if ( verbose && MyPID==0 ) {
      cout << "*** EpetraAdapter FAILED TestOperatorTraits() ***" 
           << endl << endl;
    }
    break;
  }


#ifdef TEST_THYRA_ADAPTER
  typedef Thyra::MultiVectorBase<double> TMVB;
  typedef Thyra::LinearOpBase<double>    TLOB;
  // create thyra objects from the epetra objects

  // first, a Thyra::VectorSpaceBase
  Teuchos::RefCountPtr<const Thyra::MPIVectorSpaceBase<double> > epetra_vs = 
    Thyra::create_MPIVectorSpaceBase(Map);

  // then, a ScalarProdVectorSpaceBase
  Teuchos::RefCountPtr<const Thyra::ScalarProdVectorSpaceBase<double> > sp_domain = 
    Teuchos::rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<double> >(epetra_vs,true);

  // then, a MultiVectorBase (from the Epetra_MultiVector)
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<double> > thyra_ivec = 
    Thyra::create_MPIMultiVectorBase(Teuchos::rcp_implicit_cast<Epetra_MultiVector>(ivec),epetra_vs,sp_domain);

  // then, a LinearOpBase (from the Epetra_CrsMatrix)
  Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > thyra_op = 
    Teuchos::rcp( new Thyra::EpetraLinearOp(A) );


  // test the Thyra adapter multivector
  ierr = Anasazi::TestMultiVecTraits<double,TMVB>(MyOM,thyra_ivec);
  gerr |= ierr;
  switch (ierr) {
  case Anasazi::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** ThyraAdapter PASSED TestMultiVecTraits()" << endl;
    }
    break;
  case Anasazi::Failed:
    if ( verbose && MyPID==0 ) {
      cout << "*** ThyraAdapter FAILED TestMultiVecTraits() ***" 
           << endl << endl;
    }
    break;
  }

  // test the Thyra adapter operator 
  ierr = Anasazi::TestOperatorTraits<double,TMVB,TLOB>(MyOM,thyra_ivec,thyra_op);
  gerr |= ierr;
  switch (ierr) {
  case Anasazi::Ok:
    if ( verbose && MyPID==0 ) {
      cout << "*** ThyraAdapter PASSED TestOperatorTraits()" << endl;
    }
    break;
  case Anasazi::Failed:
    if ( verbose && MyPID==0 ) {
      cout << "*** ThyraAdapter FAILED TestOperatorTraits() ***" 
           << endl << endl;
    }
    break;
  }
#endif

  // Release all objects
  delete [] NumNz;
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (gerr) {
    if (verbose && MyPID==0)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;

}
