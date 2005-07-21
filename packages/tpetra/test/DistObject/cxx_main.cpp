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
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_Util.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

template <typename OrdinalType, typename ScalarType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages);

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
	// if debug is enabled, it will still output on all nodes
	verbose = (verbose && (myImageID == 0));
  
	// start the testing
	if(verbose) outputStartMessage("DistObject");
	int ierr = 0;

	// call the actual test routines
	ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("DistObject", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages) {
	std::string className = "DistObject<" + Teuchos::OrdinalTraits<OrdinalType>::name() + "," 
		+ Teuchos::ScalarTraits<ScalarType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	int ierr = 0;
	int returnierr = 0;

	// fixtures
	
	OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
	OrdinalType const negOne = zero - one;
	OrdinalType const indexBase = zero;
	OrdinalType const numGlobalElements = intToOrdinal<OrdinalType>(10); // magic number

	// Tpetra::Platform objects
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform(MPI_COMM_WORLD);
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> esPlatform(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform;
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> esPlatform;
#endif

	// ElementSpace and VectorSpace objects
	if(verbose) cout << "Creating source Distribution..." << endl;
	// source ES: nGE elements, uniform distribution
	Tpetra::ElementSpace<OrdinalType> elementspace(numGlobalElements, indexBase, esPlatform);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platform);
	if(verbose) cout << "Creating target Distribution..." << endl;
	// dest ES: locally-replicated versions of source ES/VS
	std::vector<OrdinalType> allGIDs(numGlobalElements);
	for(OrdinalType i = zero; i < numGlobalElements; i++)
		allGIDs.at(i) = i;
	Tpetra::ElementSpace<OrdinalType> es2(negOne, numGlobalElements, allGIDs, indexBase, esPlatform);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vs2(es2, platform);

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	// ...
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	// ========================================
	// distributed -> locally replicated
	// ========================================

	if(verbose) cout << "Creating and initializing source DistObject..." << endl;
	// create the vector we'll import from
	Tpetra::Vector<OrdinalType, ScalarType> v1(vectorspace);
	// we will use column 1 from the generator
	for(OrdinalType i = zero; i < numGlobalElements; i++)
		if(vectorspace.isMyGlobalIndex(i))
			v1[vectorspace.getLocalIndex(i)] = generateValue(one, i);
	if(debug)
		cout << v1;

	if(verbose) cout << "Creating and initializing target DistObject..." << endl;
	// create the locally-replicated Vector we'll import into
	Tpetra::Vector<OrdinalType, ScalarType> v2(vs2);
	if(debug)
		cout << v2;

	if(verbose) cout << "Testing doImport... ";
	if(debug) cout << endl;
	// create the Importer and do the import
	Tpetra::Import<OrdinalType> importer(elementspace, es2);
	v2.doImport(v1, importer, Tpetra::Insert);
	if(debug) {
		cout << v1;
		cout << v2;
	}

	// extract values and check them
	std::vector<ScalarType> actualValues(v2.getNumMyEntries());
	v2.extractCopy(&actualValues.front());
	std::vector<ScalarType> expectedValues;
	generateColumn(expectedValues, one, numGlobalElements);
	/*if(debug) {
		outputData(myImageID, numImages, "Expected: " + Tpetra::toString(expectedValues));
		outputData(myImageID, numImages, "Actual: " + Tpetra::toString(actualValues));
	}*/
	if(debug && verbose)
		cout << "doImport test ";
	if(actualValues != expectedValues) {
		if(verbose)
			cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose)
			cout << "passed" << endl;
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
