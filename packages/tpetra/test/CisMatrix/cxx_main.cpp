// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

// Tpetra CisMatrix tester

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CombineMode.hpp"
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_Version.hpp"

// function prototype
template <typename OrdinalType, typename ScalarType>
int codeCoverage(bool verbose, bool debug);
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug);

int main(int argc, char* argv[]) {
	// initialize verbose & debug flags
	bool verbose = false;
	bool debug = false;
	if(argc > 1) {
		if(argv[1][0] == '-' && argv[1][1] == 'v')
			verbose = true;
		if(argv[1][0] == '-' && argv[1][1] == 'd') {
			debug = true;
			verbose = true;
		}
	}

	if (verbose)
		cout << Tpetra::Tpetra_Version() << endl << endl;
	// call test routine
	int ierr = 0;
	if(verbose) cout << "Starting CisMatrixTest..." << endl;
	//ierr += codeCoverage<int, double>(verbose, debug);
	ierr += unitTests<int, double>(verbose, debug);

	// finish up
	if(verbose)
		if(ierr == 0)
			cout << "CisMatrix test successful." << endl;
		else
			cout << "CisMatrix test failed." << endl;
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int codeCoverage(bool verbose, bool debug) {
	char const * OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
	char const * STName = Teuchos::ScalarTraits<ScalarType>::name();

	if(verbose) cout << "Starting code coverage for CisMatrix<" << OTName << "," << STName << ">." << endl;

	if(verbose) cout << "Constructors..." << endl;
	// have to create ElementSpace and VectorSpace first
  const Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;
	const Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
	Tpetra::ElementSpace<OrdinalType> elementspace(10, 0, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
  // constructor taking one VectorSpace
	Tpetra::CisMatrix<OrdinalType, ScalarType> sm(vectorspace);
  // constructor taking two VectorSpaces
  Tpetra::CisMatrix<OrdinalType, ScalarType> sm2(vectorspace, vectorspace);
	// copy constructor
	Tpetra::CisMatrix<OrdinalType, ScalarType> smClone(sm);

  // submissions
  if(verbose) cout << "Submitting entries..." << endl;
  sm.submitEntry(Tpetra::Replace, 0, 5, 1);
  sm.submitEntry(Tpetra::Replace, 0, 2, 2);
  sm.submitEntry(Tpetra::Replace, 1, 8, 0);
  sm.submitEntry(Tpetra::Replace, 1, 6, 3);
  sm.submitEntry(Tpetra::Replace, 2, 3, 2);
  sm.submitEntry(Tpetra::Replace, 3, 4, 0);
  sm.submitEntry(Tpetra::Replace, 3, 11, 1);
  sm.submitEntry(Tpetra::Replace, 3, 1, 2);
  sm.submitEntry(Tpetra::Add, 3, 1, 1);

  if(debug) cout << sm << endl;

  // scale
  if(verbose) cout << "scale..." << endl;
  sm.scale(2.0);
  if(debug) cout << sm << endl;

  // setAllToScalar
  if(verbose) cout << "setAllToScalar..." << endl;
  sm.setAllToScalar(6.0);
  if(debug) cout << sm << endl;

  // get attributes
  if(verbose) cout << "getNumGlobalNonzeros..." << endl;
  sm.getNumGlobalNonzeros(); // throw away output
  if(verbose) cout << "getNumMyNonzeros..." << endl;
  sm.getNumMyNonzeros(); // throw away output
  if(verbose) cout << "getNumGlobalDiagonals..." << endl;
  sm.getNumGlobalDiagonals(); // throw away output
  if(verbose) cout << "getNumMyDiagonals..." << endl;
  sm.getNumMyDiagonals(); // throw away output
  
  // retrieve VectorSpaces
  if(verbose) cout << "Retrieving VectorSpaces..." << endl;
  if(verbose) cout << "sm.getPrimaryDist()" << endl;
  sm.getPrimaryDist(); // throw away output
  //if(verbose) cout << "sm.getSecondaryDist()" << endl;
  //sm.getSecondaryDist(); // throw away output
  //if(verbose) cout << "sm.getDomainMap()" << endl;
  //sm.getDomainDist(); // throw away output
  //if(verbose) cout << "sm.getRangeMap()" << endl;
  //sm.getRangeDist(); // throw away output

  // fillComplete
  if(verbose) cout << "fillComplete..." << endl;
  sm.fillComplete();

  // print
  cout << sm << endl;
  

	if(verbose) cout << "Code coverage <" << OTName << ", " << STName << "> section finished." << endl;

	return(0);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
	//int ierr = 0;
	int returnierr = 0;

	if(verbose) cout << "Starting actual testing section..." << endl;
  
  // do a simple matrix-vector multiplication
  // the CisMatrix will be 4x4, and the vectors will be length 4

  // create platform/es/vs we use
  if(verbose) cout << "Creating and Initializing platform/es/ves..." << endl;
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformO;
  Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
  Tpetra::ElementSpace<OrdinalType> elementspace(4, 0, platformO);
  Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
  if(debug) {
    cout << "Output of elementspace:" << endl;
    cout << elementspace;
    cout << "Output of vectorspace:" << endl;
    cout << vectorspace;
  }
  
  // create x vector and initialize values
  if(verbose) cout << "Creating and initializing x and y vectors..." << endl;
  Tpetra::Vector<OrdinalType, ScalarType> x(vectorspace);
  x[0] = 3.0;
  x[1] = 1.0;
  x[2] = 6.0;
  x[3] = 4.0;
  // create y vector (don't need to initialize values)
  Tpetra::Vector<OrdinalType, ScalarType> y(vectorspace);
  if(debug) {
    cout << "Output of x:" << endl;
    cout << x;
    cout << "Output of y:" << endl;
    cout << y;
  }

  // create A CisMatrix and initalize values
  if(verbose) cout << "Creating A matrix..." << endl;
  Tpetra::CisMatrix<OrdinalType, ScalarType> A(vectorspace, vectorspace);
  if(debug) cout << A;
  if(verbose) cout << "Submitting values..." << endl;
  A.submitEntry(Tpetra::Insert, 0, 2.0, 0); // CombineMode, Row/Col Number, Value, Index
  A.submitEntry(Tpetra::Insert, 0, 1.0, 2); 
  A.submitEntry(Tpetra::Insert, 1, 4.0, 1); // Matrix layout is:
  A.submitEntry(Tpetra::Insert, 1, 2.0, 3); // 2 0 1 0
  A.submitEntry(Tpetra::Insert, 2, 3.0, 0); // 0 4 0 2
  A.submitEntry(Tpetra::Insert, 2, 6.0, 2); // 3 0 6 0
  A.submitEntry(Tpetra::Insert, 3, 5.0, 1); // 0 5 0 8
  A.submitEntry(Tpetra::Insert, 3, 8.0, 3);
  if(debug) cout << A;
  if(verbose) cout << "Calling fillComplete..." << endl;
  A.fillComplete();
  if(debug) cout << A;

  // output current values
  if(verbose) cout << "Finished creating & initializing." << endl;

  // call apply
  if(verbose) cout << "Calling apply..." << endl;
  A.apply(x, y);
  if(debug) {
    cout << "Output of A:" << endl;
    cout << A;
    cout << "Output of x:" << endl;
    cout << x;
    cout << "Output of y:" << endl;
    cout << y;
  }

	// finish up
	if(verbose)
		if(returnierr == 0)
			cout << "Unit tests <" 
			     << Teuchos::OrdinalTraits<OrdinalType>::name()  << ", " 
			     << Teuchos::ScalarTraits<ScalarType>::name() << "> passed." << endl;
		else
			cout << "Unit Tests <" 
			     << Teuchos::OrdinalTraits<OrdinalType>::name() << ", " 
			     << Teuchos::ScalarTraits<ScalarType>::name() << "> failed." << endl;
  
	return(returnierr);
}
