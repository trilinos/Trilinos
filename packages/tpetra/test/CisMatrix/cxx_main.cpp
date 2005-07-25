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

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_Version.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI
#include "../tpetra_test_util.hpp"

// function prototype
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

int main(int argc, char* argv[]) {
	int myImageID = 0; // assume we are on serial
	int numImages = 1; // if MPI, will be reset later
  
	// initialize MPI if needed
#ifdef TPETRA_MPI
	myImageID = -1;
	numImages = -1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numImages);
	MPI_Comm_rank(MPI_COMM_WORLD, &myImageID);
#endif // TPETRA_MPI

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

	// change verbose to only be true on Image 0
	verbose = (verbose && (myImageID == 0));
  
	// start the testing
	if(verbose) outputStartMessage("CisMatrix");
	int ierr = 0;
  
	// call the actual test routines
	ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);

	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif // TPETRA_MPI
	if(verbose) outputEndMessage("CisMatrix", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	std::string className = "CisMatrix<" + Teuchos::OrdinalTraits<OrdinalType>::name() 
		+ "," + Teuchos::ScalarTraits<ScalarType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);
	int ierr = 0;
	int returnierr = 0;

	OrdinalType const zero = intToOrdinal<OrdinalType>(0);
	OrdinalType const one = intToOrdinal<OrdinalType>(1);
	OrdinalType const two = intToOrdinal<OrdinalType>(2);
	OrdinalType const three = intToOrdinal<OrdinalType>(3);
	OrdinalType const four = intToOrdinal<OrdinalType>(4);
  
	// do a simple matrix-vector multiplication
	// the CisMatrix will be 4x4, and the vectors will be length 4

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	// ========================================
	// create platform/es/vs we use
	// ========================================
	if(verbose) cout << "Creating and Initializing platform/es/ves..." << endl;
#ifdef TPETRA_MPI
	const Tpetra::MpiPlatform<OrdinalType, OrdinalType> platformO(MPI_COMM_WORLD);
	const Tpetra::MpiPlatform<OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
	const Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformO;
	const Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
#endif // TPETRA_MPI
	Tpetra::ElementSpace<OrdinalType> elementspace(four, zero, platformO);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
	if(debug) {
		cout << "Output of elementspace:" << endl;
		cout << elementspace;
		cout << "Output of vectorspace:" << endl;
		cout << vectorspace;
	}
	
	// ========================================
	// create x vector and initialize values
	// ========================================
	// layout is: { 3 1 6 4 }
	if(verbose) cout << "Creating and initializing x and y vectors..." << endl;
	Tpetra::Vector<OrdinalType, ScalarType> x(vectorspace);
	if(vectorspace.isMyGlobalIndex(zero))
		x[vectorspace.getLocalIndex(zero)] = intToScalar<ScalarType>(3);
	if(vectorspace.isMyGlobalIndex(one))
		x[vectorspace.getLocalIndex(one)] = intToScalar<ScalarType>(1);
	if(vectorspace.isMyGlobalIndex(two))
		x[vectorspace.getLocalIndex(two)] = intToScalar<ScalarType>(6);
	if(vectorspace.isMyGlobalIndex(three))
		x[vectorspace.getLocalIndex(three)] = intToScalar<ScalarType>(4);

	// ========================================
	// create y vector (don't need to initialize values)
	// ========================================
	Tpetra::Vector<OrdinalType, ScalarType> y(vectorspace);
	if(debug) {
		if(myImageID == 0) 
			cout << "Output of x:" << endl;
		cout << x;
		elementspace.comm().barrier();
		if(myImageID == 0) 
			cout << "\nOutput of y:" << endl;
		cout << y;
		elementspace.comm().barrier();
	}

	// ========================================
	// create A CisMatrix and initalize values
	// ========================================
	// Matrix layout is:
	//  2  0  1  0
	//  0  4  0  2
	//  3  0  6  0
	//  0  5  0  8
	if(verbose) cout << "Creating A matrix..." << endl;
	Tpetra::CisMatrix<OrdinalType, ScalarType> A(vectorspace);
	if(debug) cout << A;
	if(verbose) cout << "Submitting values..." << endl;
	if(vectorspace.isMyGlobalIndex(zero)) {
		A.submitEntry(Tpetra::Insert, zero, intToScalar<ScalarType>(2), zero); // CombineMode, Row/Col Number, Value, Index
		A.submitEntry(Tpetra::Insert, zero, intToScalar<ScalarType>(1), two); 
	}
	if(vectorspace.isMyGlobalIndex(one)) {
		A.submitEntry(Tpetra::Insert, one, intToScalar<ScalarType>(4), one);
		A.submitEntry(Tpetra::Insert, one, intToScalar<ScalarType>(2), three);
	}
	if(vectorspace.isMyGlobalIndex(two)) {
		A.submitEntry(Tpetra::Insert, two, intToScalar<ScalarType>(3), zero);
		A.submitEntry(Tpetra::Insert, two, intToScalar<ScalarType>(6), two);
	}
	if(vectorspace.isMyGlobalIndex(three)) {
		A.submitEntry(Tpetra::Insert, three, intToScalar<ScalarType>(5), one);
		A.submitEntry(Tpetra::Insert, three, intToScalar<ScalarType>(8), three);
	}
	if(debug) cout << A;
	if(verbose) cout << "Calling fillComplete..." << endl;
	A.fillComplete();
	if(debug) cout << A;

	// output current values
	if(verbose) cout << "Finished creating & initializing." << endl;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	// ========================================
	// test apply
	// ========================================
	if(verbose) cout << "Testing apply... ";
	if(verbose && debug) cout << endl;
	A.apply(x, y);
	if(debug) {
		if(myImageID == 0)
			cout << "\nOutput of A:" << endl;
		elementspace.comm().barrier();
		cout << A;
		if(myImageID == 0)
			cout << "\nOutput of x:" << endl;
		elementspace.comm().barrier();
		cout << x;
		if(myImageID == 0)
			cout << "\nOutput of y:" << endl;
		elementspace.comm().barrier();
		cout << y;
	}

	// layout of y should be: { 12 12 45 37 }
	if(vectorspace.isMyGlobalIndex(zero)) {
		if(y[vectorspace.getLocalIndex(zero)] != intToScalar<ScalarType>(12))
			ierr++;
	}
	if(vectorspace.isMyGlobalIndex(one)) {
		if(y[vectorspace.getLocalIndex(one)] != intToScalar<ScalarType>(12))
			ierr++;
	}
	if(vectorspace.isMyGlobalIndex(two)) {
		if(y[vectorspace.getLocalIndex(two)] != intToScalar<ScalarType>(45))
			ierr++;
	}
	if(vectorspace.isMyGlobalIndex(three)) {
		if(y[vectorspace.getLocalIndex(three)] != intToScalar<ScalarType>(37))
			ierr++;
	}
	
	if(verbose && debug) cout << "apply test ";
	if(ierr != 0) {
		if(verbose) cout << "failed" << endl;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;

	// ========================================
	// finish up
	// ========================================
	elementspace.comm().barrier();
	if(verbose) {
		if(returnierr == 0)
			outputHeading("Unit tests for " + className + " passed.");
		else
			outputHeading("Unit tests for " + className + " failed.");
	}
	return(returnierr);
}
