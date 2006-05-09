// @HEADER
// ***********************************************************************
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
// ***********************************************************************
// @HEADER
//
//  This example computes the specified eigenvalues of the discretized 2D Laplacian
//  using the block Davidson.  

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBlockDavidson.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicEigenproblem.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  int i, info;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Anasazi::ReturnType returnCode = Anasazi::Ok; 
  int MyPID = Comm.MyPID();

  // Create an output manager to handle the I/O from the solver
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM =
    Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );
  MyOM->SetVerbosity( Anasazi::FinalSummary + Anasazi::TimingDetails );  

  std::string which;
  if (argc > 1) {
    which = argv[1];
  }
  else {
    which = "SM";
  }
  if ( which != "SM" && which != "LM" && which != "SR" && which != "LR" ) {
    if (MyOM->doPrint()) {
      std::cout << "Usage: " << argv[0] << " [sort string]" << endl
        << "where:" << endl
        << "sort string       - SM | LM | SR | LR" << endl << endl;
    }
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  //  Dimension of the matrix
  int nx = 10;                    // Discretization points in any one direction.
  int NumGlobalElements = nx*nx;  // Size of matrix nx*nx

  // Construct a Map that puts approximately the same number of
  // equations on each processor.

  Epetra_Map Map(NumGlobalElements, 0, Comm);

  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor
  std::vector<int> NumNz(NumMyElements);

  /* We are building a matrix of block structure:
  
      | T -I          |
      |-I  T -I       |
      |   -I  T       |
      |        ...  -I|
      |           -I T|

   where each block is dimension nx by nx and the matrix is on the order of
   nx*nx.  The block T is a tridiagonal matrix. 
  */

  for (i=0; i<NumMyElements; i++) {
    if (MyGlobalElements[i] == 0 || MyGlobalElements[i] == NumGlobalElements-1 || 
        MyGlobalElements[i] == nx-1 || MyGlobalElements[i] == nx*(nx-1) ) {
      NumNz[i] = 3;
    }
    else if (MyGlobalElements[i] < nx || MyGlobalElements[i] > nx*(nx-1) || 
             MyGlobalElements[i]%nx == 0 || (MyGlobalElements[i]+1)%nx == 0) {
      NumNz[i] = 4;
    }
    else {
      NumNz[i] = 5;
    }
  }

  // Create an Epetra_Matrix

  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]) );

  // Compute coefficients for discrete convection-diffution operator
  double one = 1.0;
  std::vector<double> Values(4);
  std::vector<int> Indices(4);
  double rho = 0.0;  
  double h = one /(nx+1);
  double h2 = h*h;
  double c = 5.0e-01*rho/ h;
  Values[0] = -one/h2 - c; Values[1] = -one/h2 + c; Values[2] = -one/h2; Values[3]= -one/h2;
  double diag = 4.0 / h2;
  int NumEntries;

  for (i=0; i<NumMyElements; i++)
  {
    if (MyGlobalElements[i]==0)
    {
      Indices[0] = 1;
      Indices[1] = nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == nx*(nx-1))
    {
      Indices[0] = nx*(nx-1)+1;
      Indices[1] = nx*(nx-2);
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == nx-1)
    {
      Indices[0] = nx-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
      Indices[0] = 2*nx-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
    {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
      Indices[0] = nx*(nx-1)-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] < nx)
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] > nx*(nx-1))
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i]%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]+1;
      Indices[1] = MyGlobalElements[i]-nx;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if ((MyGlobalElements[i]+1)%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]-nx;
      Indices[1] = MyGlobalElements[i]+nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
      Indices[0] = MyGlobalElements[i]-1;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      Indices[3] = MyGlobalElements[i]+nx;
      NumEntries = 4;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    // Put in the diagonal entry
    info = A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i]);
    assert( info==0 );
  }
  
  // Finish up
  info = A->FillComplete();
  assert( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
  //
  //************************************
  // Start the block Davidson iteration
  //***********************************
  //
  //  Variables used for the Block Davidson Method
  //
  int nev = 4;
  int blockSize = 5;
  int maxBlocks = 8;
  int maxIters = 300;
  double tol = 1e-8;  
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Max Iters", maxIters );
  MyPL.set( "Tol", tol );
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;
  
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();
  
  // Create the eigenproblem.
  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );
  
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true); 
  
  // Set the number of eigenvalues requested
  MyProblem->SetNEV( nev );
  
  // Inform the eigenproblem that you are finishing passing it information
  info = MyProblem->SetProblem();
  if (info) {
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;
  }

  // Create a sort manager
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySM = 
    Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );

  // Initialize the Block Davidson solver
  Anasazi::BlockDavidson<double, MV, OP> MySolver(MyProblem, MySM, MyOM, MyPL);

  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();

  // If the solver returns a code other than OK, exit the example.
  if (returnCode != Anasazi::Ok)
    return -1;

  // Retrieve eigenvalues
  Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();

  // Retrieve eigenvectors
  Teuchos::RefCountPtr<Epetra_MultiVector> evecs = MyProblem->GetEvecs();

  // Compute residuals.
  Teuchos::SerialDenseMatrix<int,double> T(nev, nev);
  Epetra_MultiVector tempAevec( Map, nev );
  std::vector<double> normA(nev);
  if (MyOM->doPrint()) {
    cout << endl;
    cout << "Actual Residuals" << endl;
    cout << "------------------------------------------------------" << endl;
  }
  T.putScalar(0.0); 
  for (i=0; i<nev; i++) {
    T(i,i) = (*evals)[i]; 
  }
  A->Apply( *evecs, tempAevec );
  MVT::MvTimesMatAddMv( -1.0, *evecs, T, 1.0, tempAevec );
  MVT::MvNorm( tempAevec, &normA );
  for (i=0; i<nev; i++)
    normA[i] /= Teuchos::ScalarTraits<double>::magnitude((*evals)[i]);

  if (MyOM->doPrint()) {
    cout.setf( ios_base::right, ios_base::adjustfield );
    cout << std::setw( 16 ) << "Eigenvalue"
         << std::setw( 16 ) << "Direct Residual" 
	 << endl;
    cout << "------------------------------------------------------" << endl;
    for (i=0; i<nev; i++) {
      cout << std::setw( 16 ) << (*evals)[i] 
	   << std::setw( 16 ) << normA[i] 
	   << endl;
    }  
    cout << "------------------------------------------------------" << endl;
  }

#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    //
    // Default return value
    //
  return 0;
}
