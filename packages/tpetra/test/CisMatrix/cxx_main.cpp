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

#include "../tpetra_test_util.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CombineMode.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI


// function prototype
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);
template <typename OrdinalType, typename ScalarType>
int testApply(bool verbose, bool debug, int myImageID, int numImages, bool isRowOriented, bool doTranspose);

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

	//mpiBreakpoint(myImageID);
  
	// call the actual test routines
	ierr += unitTests<int, float>(verbose, debug, myImageID, numImages);

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

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	// ...

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================
	
	// ========================================
	// test apply - row-oriented
	// ========================================
	if(verbose) cout << "Testing apply (row-oriented, non-transpose)... ";
	if(verbose && debug) cout << endl;
	ierr = testApply<OrdinalType, ScalarType>((verbose && debug), debug, myImageID, numImages, true, false);
	
	if(verbose && debug) cout << "apply test ";
	if(ierr != 0) {
		if(verbose) cout << "failed" << endl;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// ========================================
	// test apply - row-oriented, transpose
	// ========================================
	if(verbose) cout << "Testing apply (row-oriented, transpose)... ";
	if(verbose && debug) cout << endl;
	ierr = testApply<OrdinalType, ScalarType>((verbose && debug), debug, myImageID, numImages, true, true);
	
	if(verbose && debug) cout << "apply test ";
	if(ierr != 0) {
		if(verbose) cout << "failed" << endl;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// ========================================
	// test apply - column-oriented
	// ========================================
	if(verbose) cout << "Testing apply (column-oriented, non-transpose)... ";
	if(verbose && debug) cout << endl;
	ierr = testApply<OrdinalType, ScalarType>((verbose && debug), debug, myImageID, numImages, false, false);
	
	if(verbose && debug) cout << "apply test ";
	if(ierr != 0) {
		if(verbose) cout << "failed" << endl;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// ========================================
	// test apply - column-oriented, transpose
	// ========================================
	if(verbose) cout << "Testing apply (column-oriented, transpose)... ";
	if(verbose && debug) cout << endl;
	ierr = testApply<OrdinalType, ScalarType>((verbose && debug), debug, myImageID, numImages, false, true);
	
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
	if(verbose) {
		if(returnierr == 0)
			outputHeading("Unit tests for " + className + " passed.");
		else
			outputHeading("Unit tests for " + className + " failed.");
	}
	return(returnierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int testApply(bool verbose, bool debug, int myImageID, int numImages, bool isRowOriented, bool doTranspose) {
	OrdinalType const zero = intToOrdinal<OrdinalType>(0);
	OrdinalType const one = intToOrdinal<OrdinalType>(1);
	OrdinalType const two = intToOrdinal<OrdinalType>(2);
	OrdinalType const three = intToOrdinal<OrdinalType>(3);

	int ierr = 0;

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
	Tpetra::ElementSpace<OrdinalType> rowES(two, zero, platformO);
	Tpetra::ElementSpace<OrdinalType> colES(three, zero, platformO);
	Tpetra::VectorSpace<OrdinalType, ScalarType> primary(rowES, platformV);
	Tpetra::VectorSpace<OrdinalType, ScalarType> secondary(colES, platformV);

	// ========================================
	// create CisMatrix and initalize values
	// ========================================
	// Matrix layout is:
	//  3  2  4
	//  0  1  0
	if(verbose) cout << "Creating A matrix..." << endl;
	Teuchos::RefCountPtr< Tpetra::CisMatrix<OrdinalType, ScalarType> > A;
	if(isRowOriented)
		A = Teuchos::rcp(new Tpetra::CisMatrix<OrdinalType, ScalarType>(primary));
	else
		A = Teuchos::rcp(new Tpetra::CisMatrix<OrdinalType, ScalarType>(secondary, false));

	if(verbose) cout << "Submitting values..." << endl;
	if(isRowOriented) {
		if(primary.isMyGlobalIndex(zero)) {
			A->submitEntry(Tpetra::Insert, zero, intToScalar<ScalarType>(3), zero); // CombineMode, Row Number, Value, Col Index
			A->submitEntry(Tpetra::Insert, zero, intToScalar<ScalarType>(2), one); 
			A->submitEntry(Tpetra::Insert, zero, intToScalar<ScalarType>(4), two); 
		}
		if(primary.isMyGlobalIndex(one)) {
			A->submitEntry(Tpetra::Insert, one, intToScalar<ScalarType>(1), one);
		}
	}
	else {
		if(secondary.isMyGlobalIndex(zero)) {
			A->submitEntry(Tpetra::Insert, zero, intToScalar<ScalarType>(3), zero); // CombineMode, Col Number, Value, Row Index
		}
		if(secondary.isMyGlobalIndex(one)) {
			A->submitEntry(Tpetra::Insert, one, intToScalar<ScalarType>(2), zero); 
			A->submitEntry(Tpetra::Insert, one, intToScalar<ScalarType>(1), one); 
		}
		if(secondary.isMyGlobalIndex(two)) {
			A->submitEntry(Tpetra::Insert, two, intToScalar<ScalarType>(4), zero);
		}
	}
	
	if(verbose) cout << "Calling fillComplete..." << endl;
	A->fillComplete(secondary, primary);
	if(debug)
		cout << *A;

	// ========================================
	// create x vector and initialize values
	// ========================================
	// layout is: { 5 1 2 } for non-transpose, { 5 6 } for transpose
	if(verbose) cout << "Creating and initializing x vector..." << endl;
	Teuchos::RefCountPtr< Tpetra::Vector<OrdinalType, ScalarType> > x;
	if(!doTranspose) { // non-transpose
		x = Teuchos::rcp(new Tpetra::Vector<OrdinalType, ScalarType>(secondary));
		if(secondary.isMyGlobalIndex(zero))
			(*x)[secondary.getLocalIndex(zero)] = intToScalar<ScalarType>(5);
		if(secondary.isMyGlobalIndex(one))
			(*x)[secondary.getLocalIndex(one)] = intToScalar<ScalarType>(1);
		if(secondary.isMyGlobalIndex(two))
			(*x)[secondary.getLocalIndex(two)] = intToScalar<ScalarType>(2);
	}
	else { // transpose
		x = Teuchos::rcp(new Tpetra::Vector<OrdinalType, ScalarType>(primary));
		if(primary.isMyGlobalIndex(zero))
			(*x)[primary.getLocalIndex(zero)] = intToScalar<ScalarType>(5);
		if(primary.isMyGlobalIndex(one))
			(*x)[primary.getLocalIndex(one)] = intToScalar<ScalarType>(6);
	}
	if(debug) {
		if(verbose) 
			cout << "Output of x:" << endl;
		cout << *x;
		rowES.comm().barrier();
	}

	
	// ========================================
	// create y vector (don't need to initialize values)
	// ========================================
	// layout is: { 0 0 } for non-transpose, { 0 0 0 } for transpose
	// (y should be set to all zeros by default)
	if(verbose) cout << "Creating and initializing y vector..." << endl;
	Teuchos::RefCountPtr< Tpetra::Vector<OrdinalType, ScalarType> > y;
	if(!doTranspose)
		y = Teuchos::rcp(new Tpetra::Vector<OrdinalType, ScalarType>(primary));
	else
		y = Teuchos::rcp(new Tpetra::Vector<OrdinalType, ScalarType>(secondary));
	
	// ========================================
	// call apply
	// ========================================

	if(verbose) cout << "Calling apply..." << endl;
	A->apply(*x, *y, doTranspose);
	if(debug) {
		if(verbose) 
			cout << "\nOutput of y:" << endl;
		cout << *y;
		rowES.comm().barrier();
	}

	// ========================================
	// check results
	// ========================================
	
	if(!doTranspose) {
		// layout of y should be: { 25 1 }
		if(primary.isMyGlobalIndex(zero)) {
			if((*y)[primary.getLocalIndex(zero)] != intToScalar<ScalarType>(25))
				ierr++;
		}
		if(primary.isMyGlobalIndex(one)) {
			if((*y)[primary.getLocalIndex(one)] != intToScalar<ScalarType>(1))
				ierr++;
		}
	}
	else {
		// layout of y should be: { 15 16 20 }
		if(secondary.isMyGlobalIndex(zero)) {
			if((*y)[secondary.getLocalIndex(zero)] != intToScalar<ScalarType>(15))
				ierr++;
		}
		if(secondary.isMyGlobalIndex(one)) {
			if((*y)[secondary.getLocalIndex(one)] != intToScalar<ScalarType>(16))
				ierr++;
		}
		if(secondary.isMyGlobalIndex(two)) {
			if((*y)[secondary.getLocalIndex(two)] != intToScalar<ScalarType>(20))
				ierr++;
		}
	}

	return(ierr);
}
