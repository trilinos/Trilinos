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
#include "AnasaziModalSolverUtils.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziBasicOutputManager.hpp"

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
  double zero = 0.0;
  
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

  //  Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > om = Teuchos::rcp( new Anasazi::BasicOutputManager<double>() );

  //  Create ModalSolverUtils object
  Anasazi::ModalSolverUtils<double, Epetra_MultiVector, Epetra_Operator> msUtils( om );
  
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
  //  test sortScalars methods
  //--------------------------------------------------------------------------
  
  // Set output format
  cout.precision(6);
  cout.setf(ios::fixed, ios::floatfield);

  if (verbose && MyPID == 0)
    cout<< endl <<"*************** CHECKING sortScalars ***************" << endl << endl;
  
  // Create a vector with random entries.
  std::vector<double> Avec( NumColumns );
  std::vector<int> perm( NumColumns );
  
  if (verbose && MyPID == 0)
    cout << "***** Before sorting *****" <<endl;
  
  for (i=0; i<NumColumns; i++) {
    Avec[i] = Teuchos::ScalarTraits<double>::random();
    if (verbose && MyPID == 0)
      cout << Avec[i] << "\t";
  }	
  if (verbose && MyPID == 0)
    cout << endl << endl;
  
  std::vector<double> Avec2( Avec );
  
  if (verbose && MyPID == 0)
    cout<< endl <<"*********** sortScalars Test ***********" << endl << endl;

  // Sort vector with random entries.
  info = msUtils.sortScalars( NumColumns, &Avec[0], &perm[0] );
  if (info != 0)
    numberFailedTests++;
  
  if (verbose && MyPID == 0) {
    cout << "***** After sorting *****" <<endl;
    for (i=0; i<NumColumns; i++) {
      cout << Avec[i] << "\t";
    }	
    cout << endl;
  }
  
  // Check that the scalars are in ascending order
  double tmp = Avec[0];
  for (i=1; i<NumColumns; i++) {
    if ( Avec[i] < tmp ) {
      if (verbose && MyPID == 0)
	cout<< "ERROR:  SORTING FAILED!"<<endl;
      numberFailedTests++;
      break;
    }
  }
  
  if (verbose && MyPID == 0)
    cout<< endl <<"*********** CHECKING sortScalars_Vectors ***********" << endl << endl;
  
  Epetra_MultiVector Bvec(Map, NumColumns);
  
  Epetra_Vector* colI = 0;
  for (i=0; i<NumColumns; i++) {
    colI = Bvec(i);
    colI->PutScalar( i );
  }

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** sortScalars_Vectors Test ***********" << endl << endl;
  
  // Sort vector with random entries and corresponding multi_vector.
  info = 0;
  info = msUtils.sortScalars_Vectors( NumColumns, &Avec2[0], &Bvec );
  if (info != 0)
    numberFailedTests++;
  
  if (verbose && MyPID == 0) {
    cout << "***** After sorting *****" <<endl;
    for (i=0; i<NumColumns; i++) {
      cout << Avec2[i] << "\t";
    }	
    cout << endl << endl;
  }
  
  // Check that the scalars are in ascending order
  tmp = Avec2[0];
  for (i=1; i<NumColumns; i++) {
    if ( Avec2[i] < tmp ) {
       if (verbose && MyPID == 0)
	cout<< "ERROR:  SORTING FAILED!"<<endl;
      numberFailedTests++;
      break;
    }
  }
  
  // Check that the vectors are sorted appropriately 
  // NOTE:  Reusing Avec2
  Bvec.Norm1( &Avec2[0] );
  for (i=0; i<NumColumns; i++) {
    if ( Avec2[i]/NumGlobalElements != perm[i] ) {
      if (verbose)
	cout << Bvec << endl;;
      if (verbose && MyPID == 0)
	cout<< "ERROR:  SORTING FAILED!"<<endl;
      numberFailedTests++;
      break;
    }
  }

  
  //--------------------------------------------------------------------------
  //  test projection methods
  //--------------------------------------------------------------------------

  // Set output format
  cout.precision(6);
  cout.setf(ios::scientific, ios::floatfield);

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** CHECKING massOrthonormalize ***********" << endl << endl;

  double err = zero;

  // randomize the multivector
  Bvec.Random();

  // create mass matrix
  Epetra_CrsMatrix MM(Copy, Map, 1);
  std::vector<double> diag(NumMyElements);
  for (i=0; i<NumMyElements; i++) {
    diag[i] = Teuchos::ScalarTraits<double>::random() + 2;
    MM.InsertGlobalValues( MyGlobalElements[i], 1, &diag[i], &MyGlobalElements[i] );
  }
  int ret = MM.FillComplete();
  assert( ret==0 );
  MM.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

  // create a multivector that is the random matrix times the mass matrix
  Epetra_MultiVector MBvec(Map, Bvec.NumVectors());
  
  // multiply the Bvec by the mass matrix
  MM.Apply( Bvec, MBvec );

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** Orthonormalization Test ***********" << endl << endl;

  // orthonormalize Bvec
  info = 0;
  info = msUtils.massOrthonormalize( Bvec, MBvec, &MM, Bvec, NumColumns, 2 );

  // check error
  err = msUtils.errorOrthonormality( &Bvec, &MM );

  if (err > 1.0e-14 || info != 0) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec << endl;
    if (verbose && MyPID == 0) {
      cout<< "massOrthonormalize return code: "<< info << endl;
      cout<< "Orthonormality error: "<< err << endl;
      cout<< "ERROR:  ORTHONORMALITY TEST FAILED [MASS ORTHONORMALIZATION]!"<<endl;
    }
   } 
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err << endl;
      cout<< "ORTHONORMALITY TEST PASSED [MASS ORTHONORMALIZATION]!"<<endl;
    }
  }    
  
  // create and randomize some test matrices
  Epetra_MultiVector Bvec_small( Map, NumColumns - 2 ), MBvec_small( Map, NumColumns - 2 );
  Bvec_small.Random();

  // multiply Bvec_small by the mass matrix 
  MM.Apply( Bvec_small, MBvec_small );

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** Orthogonalization Test ***********" << endl << endl;

  // orthogonalize Bvec_small against Bvec
  info = 0;
  info = msUtils.massOrthonormalize( Bvec_small, MBvec_small, &MM, Bvec, NumColumns-2, 1 );

  // check error
  err = msUtils.errorOrthogonality( &Bvec_small, &Bvec, &MM );  

  if (err > 1.0e-14 || info != 0) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec_small << endl;
    if (verbose && MyPID == 0) {
      cout<< "massOrthonormalize return code: "<< info << endl;
      cout<< "Orthogonality error: "<< err << endl;
      cout<< "ERROR:  ORTHOGONALITY TEST FAILED [MASS ORTHOGONALIZATION]!"<<endl;
    }
  }
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthogonality error: "<< err << endl;
      cout<< "ORTHOGONALITY TEST PASSED [MASS ORTHOGONALIZATION]!"<<endl;
    }
  }    
  
  // check orthonormality of  Bvec_small, should be bad.
  err = msUtils.errorOrthonormality( &Bvec_small );
  if (verbose && MyPID == 0)
    cout<< "Orthonormality error: "<< err << " (should not be small) "<< endl;

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** Orthonormalization Combo Test ***********" << endl << endl;

  Epetra_MultiVector Bvec_big( Map, NumColumns + 2 );
  Bvec_big.Random();

  // randomize the multivector
  Bvec.Random();

  // orthonormalize Bvec
  info = 0;
  info = msUtils.massOrthonormalize( Bvec, Bvec, 0, Bvec, NumColumns, 2 );

  // check error
  err = msUtils.errorOrthonormality( &Bvec );

  if (err > 1.0e-14 || info != 0) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec << endl;
    if (verbose && MyPID == 0) {
      cout<< "massOrthonormalize return code: "<< info << endl;
      cout<< "Orthonormality error: "<< err << endl;
      cout<< "ERROR:  ORTHO COMBO TEST FAILED [NO MASS ORTHOGONALIZATION]!"<<endl;
    }
   } 
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err << endl;
      cout<< "ORTHO COMBO TEST PASSED [NO MASS ORTHOGONALIZATION]!"<<endl;
    }
  }    

  // orthonormalize Bvec_large against Bvec
  info = 0;
  info = msUtils.massOrthonormalize( Bvec_big, Bvec_big, 0, Bvec, NumColumns+2, 0 );

  // check error
  double err1 = msUtils.errorOrthonormality( &Bvec_big );
  double err2 = msUtils.errorOrthogonality( &Bvec_big, &Bvec );  

  if ((err1 > 1.0e-14) || (err2 > 1.0e-14) || info != 0) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec_big << endl;
    if (verbose && MyPID == 0) {
      cout<< "massOrthonormalize return code: "<< info << endl;
      cout<< "Orthonormality error: "<< err1 << endl;
      cout<< "Orthogonality error: "<< err2 << endl;
      cout<< "ERROR:  ORTHO COMBO TEST FAILED [NO MASS ORTHONORMALIZATION]!"<<endl;
    }
  }
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err1 << endl;
      cout<< "Orthogonality error: "<< err2 << endl;
      cout<< "ORTHO COMBO TEST PASSED [NO MASS ORTHONORMALIZATION]!"<<endl;
    }
  }    
  
  // check direct solver code

  if (verbose && MyPID == 0)
    cout<< endl <<"************* Direct Solver Test *************" << endl << endl;

  // using orthogonal vectors in Bvec in previous test to project the
  // matrix with sorted vectors in Avec as diagonal entries.

  int nev = ((NumColumns - 3) > 0) ? NumColumns - 3 : 3;
  std::vector<double> lambda( nev );
  Teuchos::SerialDenseMatrix<int,double> K( NumColumns, NumColumns );
  Teuchos::SerialDenseMatrix<int,double> M( NumColumns, NumColumns );
  Teuchos::SerialDenseMatrix<int,double> EV( NumColumns, nev );

  for (i=0; i<NumColumns; i++) {
    M(i,i) = 1.0;
    K(i,i) = Avec[i];
  }
  
  Anasazi::MultiVecTraits<double,MV>::MvTimesMatAddMv( 1.0, Bvec, K, 0.0, MBvec );
  Anasazi::MultiVecTraits<double,MV>::MvTransMv( 1.0, Bvec, MBvec, K );

  // Compute eigenvalues
  info = 0;
  info = msUtils.directSolver( NumColumns, K, 0, &EV, &lambda, &nev, 10 );

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
  msUtils.sortScalars( NumColumns, &true_lambda[0] );

  // Compute eigenvalues
  info = 0;
  info = msUtils.directSolver( NumColumns, K, &M, &EV, &lambda, &nev, 1 );

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
  info = msUtils.directSolver( NumColumns, K, &M, &EV, &lambda, &nev, 0 );

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
