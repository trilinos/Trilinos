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
#include "Tpetra_BasicDirectory.hpp"
#include "Tpetra_Util.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

template <typename OrdinalType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages);

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
  
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = (verbose && (myImageID == 0));
  
  // start the testing
	if(verbose) outputStartMessage("Directory");
  int ierr = 0;
  
  // call the actual test routines
  ierr += unitTests<int>(verbose, debug, myImageID, numImages);
	//ierr += unitTests<unsigned int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Directory", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages) {
  std::string className = "Directory<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;

  // fixtures
  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const indexBase = zero;
  OrdinalType numGlobalElements;
  std::vector<OrdinalType> allGIDs;
  std::vector<int> imageIDs;
  std::vector<int> expectedImageIDs;
  std::vector<OrdinalType> LIDs;
  std::vector<OrdinalType> expectedLIDs;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  
  // create platform needed for directory construction
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
#endif
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  // ========================================
  // test with a uniform ES
  // ========================================
  if(verbose) cout << "Testing Directory using a uniform contiguous ElementSpace (ES ctr 1)... " << endl;
  numGlobalElements = intToOrdinal<OrdinalType>(10); // 10 elements, distributed uniformly by ES
  Tpetra::ElementSpace<OrdinalType> elementspace(numGlobalElements, indexBase, platform);
  Tpetra::BasicDirectory<OrdinalType> directory(elementspace);

  // fill allGIDs with {0, 1, 2... 9}
  for(OrdinalType i = zero; i < numGlobalElements; i++)
    allGIDs.push_back(i);

  // fill expectedImageIDs with what should be the values
  // divide numGlobalElements evenly, give remainder to first images (one each)
  expectedImageIDs.reserve(numGlobalElements);
  int numEach = numGlobalElements / numImages;
  int remainder = numGlobalElements - (numEach * numImages);
  for(int i = 0; i < numImages; i++) {
    for(int j = 0; j < numEach; j++)
      expectedImageIDs.push_back(i);
    if(remainder > 0) {
      expectedImageIDs.push_back(i);
      remainder--;
    }
  }
  
  if(verbose) cout << "Testing getDirectoryEntries(imageIDs only)... ";
  directory.getDirectoryEntries(allGIDs, imageIDs);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "GIDs: " + Tpetra::toString(allGIDs));
    outputData(myImageID, numImages, "imageIDs: " + Tpetra::toString(imageIDs));
    outputData(myImageID, numImages, "Expected: " + Tpetra::toString(expectedImageIDs));
    if(verbose) cout << "getDirectoryEntries(imageIDs only) test ";
  }
  if(imageIDs != expectedImageIDs) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // create expectedLIDs array
  expectedLIDs.reserve(numGlobalElements);
  int curImageID = expectedImageIDs.front();
  int curLID = 0;
  for(int i = 0; i < numGlobalElements; i++) { // go through entire LIDs array
    if(expectedImageIDs[i] > curImageID) {
      curLID = 0;
      curImageID = expectedImageIDs[i];
    }
    expectedLIDs.push_back(curLID);
    curLID++;
  }

  // copy old imageID results
  std::vector<int> prevImageIDs = imageIDs;

  if(verbose) cout << "Testing getDirectoryEntries(imageIDs and LIDs)... ";
  if(imageIDs != prevImageIDs) {
    outputData(myImageID, numImages, "\nERROR: imageIDs returned by two getDirectoryEntries functions do not match");
  }
  directory.getDirectoryEntries(allGIDs, imageIDs, LIDs);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "GIDs: " + Tpetra::toString(allGIDs));
    outputData(myImageID, numImages, "LIDs: " + Tpetra::toString(LIDs));
    outputData(myImageID, numImages, "Expected: " + Tpetra::toString(expectedLIDs));
    if(verbose) cout << "getDirectoryEntries(imageIDs and LIDs) test ";
  }
  if(LIDs != expectedLIDs) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
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
