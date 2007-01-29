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
//  This test is for the internal utilities that are used by the modal analysis solver.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziSolverUtils.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) 
{
  int i, info = 0;
  
#ifdef EPETRA_MPI
  
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  Epetra_SerialComm Comm;

#endif
  
  int MyPID = Comm.MyPID();
  
  bool testFailed = false;
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  
  if (verbose && MyPID == 0)
    cout << Anasazi::Anasazi_Version() << endl << endl;
  
  int numberFailedTests = 0;

  //  Create SolverUtils object
  typedef Anasazi::SolverUtils<double, Epetra_MultiVector, Epetra_Operator> Utils;
  
  //  Dimension of the multivector
  int NumGlobalElements = 99;
  int NumColumns = 7;
  
  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(NumGlobalElements, 0, Comm);
  
  int NumMyElements = Map.NumMyElements();
  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);
  
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;

  //--------------------------------------------------------------------------
  //  test direct solver code
  //--------------------------------------------------------------------------

  if (verbose && MyPID == 0)
    cout<< endl <<"************* Direct Solver Test *************" << endl << endl;

  // using orthogonal vectors in Bvec in previous test to project the
  // matrix with sorted vectors in Avec as diagonal entries.

  std::vector<double> Avec( NumColumns );
  Epetra_MultiVector Bvec(Map, NumColumns);
  Epetra_Vector* colI = 0;
  for (i=0; i<NumColumns; i++) {
    colI = Bvec(i);
    colI->PutScalar( i );
  }
  
  for (i=0; i<NumColumns; i++) {
    Avec[i] = Teuchos::ScalarTraits<double>::random();
  }
  std::sort(Avec.begin(),Avec.end());
  if (verbose && MyPID == 0) {
    for (i=0; i<NumColumns; i++) {
      cout << Avec[i] << "\t";
    }
    cout << endl;
  }

  int nev = ((NumColumns - 3) > 0) ? NumColumns - 3 : 3;
  std::vector<double> lambda( nev );
  Teuchos::SerialDenseMatrix<int,double> K( NumColumns, NumColumns );
  Teuchos::SerialDenseMatrix<int,double> M( NumColumns, NumColumns );
  Teuchos::SerialDenseMatrix<int,double> EV( NumColumns, nev );

  for (i=0; i<NumColumns; i++) {
    M(i,i) = 1.0;
    K(i,i) = Avec[i];
  }

  Epetra_MultiVector MVtmp(Map, Bvec.NumVectors());

  // orthonormalize this vector
  Anasazi::BasicOrthoManager<double,MV,OP> OM;
  OM.normalize(Bvec,Teuchos::null,Teuchos::null);

  Anasazi::MultiVecTraits<double,MV>::MvTimesMatAddMv( 1.0, Bvec, K, 0.0, MVtmp );
  Anasazi::MultiVecTraits<double,MV>::MvTransMv( 1.0, Bvec, MVtmp, K );

  // Compute eigenvalues
  info = 0;
  info = Utils::directSolver( NumColumns, K, 0, &EV, &lambda, &nev, 10 );

  testFailed = false;
  if (info != 0) {
    numberFailedTests++;
    if (verbose && MyPID == 0) {
      cout<< "directSolver return code: "<< info << endl;
      cout<< "ERROR:  DIRECT SOLVER FAILED [esType = 10]!"<<endl;
    }
  }
  else {
    for (i=0; i<nev; i++) {
      if ( Teuchos::ScalarTraits<double>::magnitude( lambda[i]-K(i,i) ) > 1.0e-14 ) {
        numberFailedTests++;
        cout<< "ERROR (Processor "<<MyPID<<"):  DIRECT SOLVER FAILED [esType = 10]!"<<endl;
        testFailed = true;
        break;
      }      
    }
    if (!testFailed && (verbose && MyPID == 0)) {
      cout<< "DIRECT SOLVER PASSED [esType = 10]!"<<endl;
    }
  }

  // now use non identity mass matrix and see how it goes

  std::vector<double> Mdiag(NumColumns);
  std::vector<double> true_lambda(NumColumns);
  for (i=0; i<NumColumns; i++) {
    Mdiag[i] = Teuchos::ScalarTraits<double>::random() + 2.0;
    M(i,i) = Mdiag[i];
    true_lambda[i] = K(i,i) / Mdiag[i]; 
  }
  //Utils::sortScalars( NumColumns, &true_lambda[0] );

  // Compute eigenvalues
  info = 0;
  info = Utils::directSolver( NumColumns, K, &M, &EV, &lambda, &nev, 1 );

  testFailed = false;
  if (info != 0) {
    if (verbose && MyPID == 0) {
      cout<< "directSolver return code: "<< info << endl;
      cout<< "ERROR:  DIRECT SOLVER FAILED [esType = 1]!"<<endl;
    }
  }
  else {
    for (i=0; i<nev; i++) {
      if ( Teuchos::ScalarTraits<double>::magnitude( lambda[i]-true_lambda[i] ) > 1.0e-14 ) {
        numberFailedTests++;
        cout<< "ERROR (Processor "<<MyPID<<"):  DIRECT SOLVER FAILED [esType = 1]!"<<endl;
        testFailed = true;
        break;
      } 
    }
    if (!testFailed && (verbose && MyPID == 0)) {
      cout<< "DIRECT SOLVER PASSED [esType = 1]!"<<endl;
    }
  }
  
  // same test using deflated eigensolver

  // Compute eigenvalues
  info = 0;
  info = Utils::directSolver( NumColumns, K, &M, &EV, &lambda, &nev, 0 );

  testFailed = false;
  if (info != 0) {
    if (verbose && MyPID == 0) {
      cout<< "directSolver return code: "<< info << endl;
      cout<< "ERROR:  DIRECT SOLVER FAILED [esType = 0]!"<<endl;
    }
  }
  else {
    for (i=0; i<nev; i++) {
      if ( Teuchos::ScalarTraits<double>::magnitude( lambda[i]-true_lambda[i] ) > 1.0e-14 ) {
        numberFailedTests++;
        cout<< "ERROR (Processor "<<MyPID<<"):  DIRECT SOLVER FAILED [esType = 0]!"<<endl;
        testFailed = true;
        break;
      }      
    }
    if (!testFailed && (verbose && MyPID == 0)) {
      cout<< "DIRECT SOLVER PASSED [esType = 0]!"<<endl;
    }
  }
  
#ifdef EPETRA_MPI

  MPI_Finalize() ;

#endif

 if (numberFailedTests) {
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
