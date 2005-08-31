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
#include <Teuchos_Array.hpp>
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_BlockElementSpace.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

int main(int argc, char* argv[]) {
	int myImageID = 0; // assume we are on serial
	int numImages = 1; // if MPI, will be reset later
  
	// initialize MPI if needed
#ifdef TPETRA_MPI
	numImages = -1;
	myImageID = -1;
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
	if(verbose) outputStartMessage("BlockElementSpace");
	int ierr = 0;
  
	// call the actual test routines
	ierr += unitTests<int>(verbose, debug, myImageID, numImages);

	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("BlockElementSpace", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	std::string className = "BlockElementSpace<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	int ierr = 0;
	int returnierr = 0;

	OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
	OrdinalType const negOne = zero - one;

	// fixtures
	OrdinalType const numElements = intToOrdinal<OrdinalType>(5);
	OrdinalType const indexBase = intToOrdinal<OrdinalType>(0);
	OrdinalType const elementSize = intToOrdinal<OrdinalType>(2);
	bool result = false;
  
	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	// create Platform
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
#endif // TPETRA_MPI

	// create ElementSpace objects (3 of them, 1 for each ES ctr)
	Tpetra::ElementSpace<OrdinalType> es1(numElements, indexBase, platform);
	Tpetra::ElementSpace<OrdinalType> es2(negOne, numElements, indexBase, platform);
	Teuchos::Array<OrdinalType> gidList;
	generateColumn(gidList, myImageID, numElements);
	Tpetra::ElementSpace<OrdinalType> es3(negOne, numElements, gidList, indexBase, platform);

	// fixed-size ctr
	if(verbose) cout << "BlockElementSpace constructor (fixed-size)...";
	Tpetra::BlockElementSpace<OrdinalType> bes1(es1, elementSize);
	if(debug) 
		cout << bes1;

	// variable-sized ctr
	if(verbose) cout << "BlockElementSpace constructor (variable-sized)...";
	Teuchos::Array<OrdinalType> eSizeList = Teuchos::tuple(one, one, one+one, one+one+one, one+one+one+one+one);
	Tpetra::BlockElementSpace<OrdinalType> bes2(es1, eSizeList);
	if(debug) cout << bes2;

	// copy ctr
	if(verbose) cout <<"BlockElementSpace copy constructor...";
	Tpetra::BlockElementSpace<OrdinalType> bes3(bes1);
	if(debug) cout << bes3;

	if(verbose) cout << "Creating compatible ElementSpace..." << endl;
	Tpetra::ElementSpace<OrdinalType> const* bes2es = bes2.generateCompatibleElementSpace();
	if(debug) cout << (*bes2es);
	delete bes2es;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	// ========================================
	// isSameAs
	// ========================================
	if(verbose) cout << "Testing isSameAs... ";
	Tpetra::BlockElementSpace<OrdinalType> bes4(es1, elementSize + one);
	result = bes1.isSameAs(bes3);
	if(result == false)
		ierr++;
	result = bes1.isSameAs(bes4);
	if(result == true)
		ierr++;

	if(verbose) {
		if(ierr != 0) 
			cout << "Failed" << endl;
		else 
			cout << "Passed" << endl;
	}
	returnierr += ierr;
	ierr = 0;

	// ========================================
	// assignment operator
	// ========================================
	if(verbose) cout << "Checking assignment operator...";
	result = bes1.isSameAs(bes4);
	if(result == true)
		ierr++;
	bes4 = bes1;
	result = bes1.isSameAs(bes4);
	if(result == false)
		ierr++;

	if(verbose) {
		if(ierr != 0) 
			cout << "Failed" << endl;
		else 
			cout << "Passed" << endl;
	}
	returnierr += ierr;
	ierr = 0;

	// ======================================================================
	// finish up
	// ======================================================================
  
	if(verbose) {
		if(returnierr == 0)
			outputHeading("Unit tests for " + className + " passed.");
		else
			outputHeading("Unit tests for " + className + " failed.");
	}
	return(returnierr);
}
