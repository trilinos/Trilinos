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
#include "Epetra_Vector.h"

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
  
#endif
  
#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  
  bool verbose = 0;
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;
  
  //	if (verbose)
  //	  cout << Anasazi::Anasazi_Version() << endl << endl;
  
  int numberFailedTests = 0;
  int returnCode = 0;  

  //  Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > om = Teuchos::rcp( new Anasazi::OutputManager<double>() );

  //  Create ModalSolverUtils object
  Anasazi::ModalSolverUtils<double, Epetra_MultiVector, Epetra_Operator> msUtils( om );
  
  //  Dimension of the multivector
  int NumGlobalElements = 99;
  int NumColumns = 7;
  
  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  Epetra_Map Map(NumGlobalElements, 0, Comm);
  
  typedef Anasazi::MultiVec<double> MV;
  typedef Anasazi::Operator<double> OP;

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

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** Orthonormalization Test ***********" << endl << endl;

  // orthonormalize Bvec
  msUtils.massOrthonormalize( Bvec, Bvec, 0, Bvec, NumColumns, 2 );

  // check error
  err = msUtils.errorOrthonormality( &Bvec );

  if (err > 1.0e-14) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec << endl;
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err << endl;
      cout<< "ERROR:  ORTHONORMALITY TEST FAILED!"<<endl;
    }
   } 
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err << endl;
      cout<< "ORTHONORMALITY TEST PASSED!"<<endl;
    }
  }    
  
  // create and randomize some test matrices
  Epetra_MultiVector Bvec_big( Map, NumColumns + 2 );
  Bvec_big.Random();
  Epetra_MultiVector Bvec_small( Map, NumColumns - 2 );
  Bvec_small.Random();
  
  if (verbose && MyPID == 0)
    cout<< endl <<"*********** Orthogonalization Test ***********" << endl << endl;

  // orthogonalize Bvec_small against Bvec
  msUtils.massOrthonormalize( Bvec_small, Bvec_small, 0, Bvec, NumColumns-2, 1 );

  // check error
  err = msUtils.errorOrthogonality( &Bvec_small, &Bvec );  

  if (err > 1.0e-14) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec_small << endl;
    if (verbose && MyPID == 0) {
      cout<< "Orthogonality error: "<< err << endl;
      cout<< "ERROR:  ORTHOGONALITY TEST FAILED!"<<endl;
    }
  }
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthogonality error: "<< err << endl;
      cout<< "ORTHOGONALITY TEST PASSED!"<<endl;
    }
  }    
  
  // check orthonormality of  Bvec_small, should be bad.
  err = msUtils.errorOrthonormality( &Bvec_small );
  if (verbose && MyPID == 0)
    cout<< "Orthonormality error: "<< err << " (should not be small) "<< endl;

  if (verbose && MyPID == 0)
    cout<< endl <<"*********** Orthonormalization Combo Test ***********" << endl << endl;

  // orthonormalize Bvec_large against Bvec
  msUtils.massOrthonormalize( Bvec_big, Bvec_big, 0, Bvec, NumColumns+2, 0 );

  // check error
  double err1 = msUtils.errorOrthonormality( &Bvec_big );
  double err2 = msUtils.errorOrthogonality( &Bvec_big, &Bvec );  

  if ((err1 > 1.0e-14) || (err2 > 1.0e-14)) {
    numberFailedTests++;
    if (verbose)
      cout<< Bvec_big << endl;
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err1 << endl;
      cout<< "Orthogonality error: "<< err2 << endl;
      cout<< "ERROR:  ORTHO COMBO TEST FAILED!"<<endl;
    }
  }
  else {
    if (verbose && MyPID == 0) {
      cout<< "Orthonormality error: "<< err1 << endl;
      cout<< "Orthogonality error: "<< err2 << endl;
      cout<< "ORTHO COMBO TEST PASSED!"<<endl;
    }
  }    
  
  // Check to see if any tests have failed
  if (numberFailedTests)
    return numberFailedTests;
  
  return 0;
}	
