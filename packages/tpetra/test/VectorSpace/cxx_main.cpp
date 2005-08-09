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
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

template <typename OrdinalType, typename ScalarType>
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
	if(verbose) outputStartMessage("VectorSpace");
	int ierr = 0;
  
	// call the actual test routines
	ierr += unitTests<int, float>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);

	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("VectorSpace", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	std::string className = "Import<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ", " + Teuchos::ScalarTraits<ScalarType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	int ierr = 0;
	int returnierr = 0;

	OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	OrdinalType const negOne = zero - Teuchos::OrdinalTraits<OrdinalType>::one();

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
	
	// create Platform and Comm
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
	Tpetra::MpiComm<OrdinalType, ScalarType> comm(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;
	Tpetra::SerialComm<OrdinalType, ScalarType> comm;
#endif
	// create ElementSpace and BlockElementSpace
	OrdinalType const numGlobalElements = intToOrdinal<OrdinalType>(10);
	OrdinalType const indexBase = intToOrdinal<OrdinalType>(2);
	OrdinalType const numPoints = intToOrdinal<OrdinalType>(3);
	Tpetra::ElementSpace<OrdinalType> es(numGlobalElements, indexBase, platformE);
	Tpetra::BlockElementSpace<OrdinalType> bes(es, numPoints);

	// constructors
	if(verbose) cout << "VectorSpace constructor (taking an ElementSpace)..." << endl;
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(es, platformV);
 
	if(verbose) cout << "VectorSpace constructor (taking a BlockElementSpace)..." << endl;
	Tpetra::VectorSpace<OrdinalType, ScalarType> blockvectorspace(bes, platformV);

	if(verbose) cout << "VectorSpace copy constructor..." << endl;
	Tpetra::VectorSpace<OrdinalType, ScalarType> v2(vectorspace);

	// print
	if(debug) {
		if(verbose) cout << "Overloaded << operator..." << endl;
		cout << vectorspace << endl;
	}
	
	// attribute access
	if(verbose) cout << "Attribute access methods..." << endl;
	OrdinalType temp = 0;
	temp = vectorspace.getNumGlobalEntries();
	temp = vectorspace.getNumMyEntries();
	temp = vectorspace.getIndexBase();
	temp = vectorspace.getMinLocalIndex();
	temp = vectorspace.getMaxLocalIndex();
	temp = vectorspace.getMinGlobalIndex();
	temp = vectorspace.getMaxGlobalIndex();
	if(vectorspace.getNumMyEntries() > 0) {
		temp = 0;
		temp = vectorspace.getGlobalIndex(temp);
		temp = vectorspace.getLocalIndex(temp);
		bool tempB = vectorspace.isMyLocalIndex(temp);
		tempB = vectorspace.isMyGlobalIndex(temp);
	}
  
	// vector creation
	if(verbose) cout << "Vector creation..." << endl;
	Tpetra::Vector<OrdinalType, ScalarType>* vecptr = vectorspace.createVector();
	temp = vecptr->getNumMyEntries();
	delete vecptr;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	bool result = false; // fixture

	// ========================================
	// isCompatible
	// ========================================
	if(verbose) cout << "Testing isCompatible... ";
	if(verbose && debug) cout << endl;
	const Tpetra::ElementSpace<OrdinalType> compatibleES(intToOrdinal<OrdinalType>(10),
														 intToOrdinal<OrdinalType>(0),
														 platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> compatibleVS(compatibleES, platformV);
	result = vectorspace.isCompatible(compatibleVS);
	if(debug)
		outputData(myImageID, numImages, "with compatibleVS: " + Tpetra::toString(result));
	if(result == false)
		ierr++;
	const Tpetra::ElementSpace<OrdinalType> differentES(intToOrdinal<OrdinalType>(15),
														intToOrdinal<OrdinalType>(2),
														platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> differentVS(differentES, platformV);
	result = differentVS.isCompatible(vectorspace);
	if(debug)
		outputData(myImageID, numImages, "with differentVS: " + Tpetra::toString(result));
	if(result == true)
		ierr++;

	if(verbose && debug) cout << "isCompatible test ";
	if(verbose) {
		if(ierr != 0) 
			cout << "Failed" << endl;
		else 
			cout << "Passed" << endl;
	}
	returnierr += ierr;
	ierr = 0;
	
	// ========================================
	// isSameAs
	// ========================================
	if(verbose) cout << "Testing isSameAs... ";
	result = vectorspace.isSameAs(v2);
	if(result == false)
		ierr++;
	result = vectorspace.isSameAs(differentVS);
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
	if(verbose) cout << "Testing asignment operator... ";
	result = vectorspace.isSameAs(differentVS);
	if(result == true)
		ierr++;
	differentVS = vectorspace;
	result = vectorspace.isSameAs(differentVS);
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
