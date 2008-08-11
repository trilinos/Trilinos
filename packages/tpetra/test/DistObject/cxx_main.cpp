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
#include "Tpetra_Import.hpp"
#include "Tpetra_CombineMode.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI
#include <Teuchos_Array.hpp>

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

	//mpiBreakpoint(myImageID);

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

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	// Tpetra::Platform objects
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform(MPI_COMM_WORLD);
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> esPlatform(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform;
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> esPlatform;
#endif
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	// fixtures
	OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
	OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
	OrdinalType const negOne = zero - one;
	OrdinalType indexBase = zero;
	OrdinalType numGlobalElements = negOne; // dummy value
	OrdinalType numMyElements = negOne; // dummy value
	std::vector<OrdinalType> myGIDs;
	std::vector<ScalarType> expectedValues;
	std::vector<ScalarType> actualValues;
	ScalarType const scalarZero = Teuchos::ScalarTraits<ScalarType>::zero();

	// ========================================
	// Tpetra::Vector, importing
	// uniquely owned -> replicated
	// ========================================
	if(verbose) cout << "Testing Tpetra::Vector, importing... ";
	if(verbose && debug) cout << endl;
	// in source, each image will have 3 elements, uniform distribution
	// in target, each image will have all of the elements
	
	// create Source ES, VS, and Vector
	if(debug && verbose) cout << "Creating and initializing source DistObject..." << endl;
	numGlobalElements = intToOrdinal<OrdinalType>(3 * numImages);
	Tpetra::ElementSpace<OrdinalType> uniqueES(numGlobalElements, indexBase, esPlatform);
	Tpetra::VectorSpace<OrdinalType, ScalarType> uniqueVS(uniqueES, platform);
	Tpetra::Vector<OrdinalType, ScalarType> v1(uniqueVS);
	// we will use column 1 from the generator for the actual values
	for(OrdinalType i = zero; i < numGlobalElements; i++)
		if(uniqueVS.isMyGlobalIndex(i))
			v1[uniqueVS.getLocalIndex(i)] = generateValue(one, i);
	if(debug)
		cout << v1;

	// create Target ES, VS, and Vector
	if(debug && verbose) cout << "Creating and initializing target DistObject..." << endl;
	numMyElements = numGlobalElements; // nME in target = nGE in source
	myGIDs.resize(numMyElements);
	for(OrdinalType i = zero; i < numMyElements; i++)
		myGIDs[i] = i;
	Tpetra::ElementSpace<OrdinalType> replicatedES(negOne, numMyElements, myGIDs, indexBase, esPlatform);
	Tpetra::VectorSpace<OrdinalType, ScalarType> replicatedVS(replicatedES, platform);
	Tpetra::Vector<OrdinalType, ScalarType> v2(replicatedVS);
	v2.setAllToScalar(scalarZero);
	if(debug)
		cout << v2;

	if(debug && verbose) cout << "Calling doImport(using importer)... " << endl;
	// create the Importer and do the import
	Tpetra::Import<OrdinalType> importer(uniqueES, replicatedES);
	v2.doImport(v1, importer, Tpetra::Insert);

	// extract values and check them
	actualValues.resize(v2.getNumMyEntries());
	v2.extractCopy(&actualValues.front());
	generateColumn(expectedValues, one, numMyElements);
	if(debug) {
		outputData(myImageID, numImages, "Expected: " + Tpetra::toString(expectedValues));
		outputData(myImageID, numImages, "Actual: " + Tpetra::toString(actualValues));
	}
	if(debug && verbose)
		cout << "Import test ";
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

	// ========================================
	// Tpetra::Vector, exporting
	// replicated -> uniquely owned
	// ========================================
	if(verbose) cout << "Testing Tpetra::Vector, exporting... ";
	if(verbose && debug) cout << endl;
	// now we use same source and target vectors from previous test,
	// but use them to do a reverse export (using the same importer)

	// initialize v3 (our new source) to values from generator
	if(debug && verbose) cout << "Writing source values..." << endl;
	Tpetra::Vector<OrdinalType, ScalarType> v3(replicatedVS);
	for(OrdinalType i = zero; i < v3.getNumMyEntries(); i++)
		v3[i] = generateValue<ScalarType>(myImageID, i);
	if(debug)
		cout << v3;

	// initialize v4 (our new target) to all zeros
	if(debug && verbose) cout << "Clearing target values..." << endl;
	Tpetra::Vector<OrdinalType, ScalarType> v4(uniqueVS);
	v4.setAllToScalar(scalarZero);
	if(debug)
		cout << v4;

	// do the export
	if(debug && verbose) cout << "Calling doExport(using exporter)... " << endl;
	Tpetra::Export<OrdinalType> exporter(replicatedES, uniqueES);
	v4.doExport(v3, exporter, Tpetra::Add);

	// extract values and check them
	actualValues.resize(v4.getNumMyEntries());
	v4.extractCopy(&actualValues.front());
	generateRowSums(expectedValues, 0, (numImages - 1), uniqueES.getMinMyGID(), uniqueES.getMaxMyGID());
	if(debug) {
		outputData(myImageID, numImages, "Expected: " + Tpetra::toString(expectedValues));
		outputData(myImageID, numImages, "Actual: " + Tpetra::toString(actualValues));
	}
	if(debug && verbose)
		cout << "Export test ";
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
	
	// ========================================
	// Tpetra::Vector, reverse importing
	// uniquely owned -> replicated
	// ========================================
	/*
	// *** THIS TEST IS NOT FINISHED YET ***
	Tpetra::ElementSpace<OrdinalType> sourceES((one + one), zero, esPlatform);
	Tpetra::VectorSpace<OrdinalType, ScalarType> sourceVS(sourceES, platform);
	Teuchos::Array<OrdinalType> targetGIDs;
	if(myImageID == 0)
		targetGIDs = Teuchos::tuple(zero, one);
	else if(myImageID == 1)
		targetGIDs = Teuchos::tuple(zero);
	Tpetra::ElementSpace<OrdinalType> targetES(negOne, targetGIDs.size(), targetGIDs, zero, esPlatform);
	Tpetra::VectorSpace<OrdinalType, ScalarType> targetVS(targetES, platform);

	Tpetra::Export<OrdinalType> temp(targetES, sourceES);
	exporter = temp;
	
	Tpetra::Vector<OrdinalType, ScalarType> target(targetVS);
	target.setAllToScalar(scalarZero);
	Tpetra::Vector<OrdinalType, ScalarType> source(sourceVS);
	OrdinalType index = zero;
	ScalarType value;
	if(myImageID == 0)
		value = intToScalar<ScalarType>(5);
	else if(myImageID == 1)
		value = intToScalar<ScalarType>(6);
	source.submitEntries(one, &index, &value);

	target.doImport(source, exporter, Tpetra::Insert);
	*/

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
