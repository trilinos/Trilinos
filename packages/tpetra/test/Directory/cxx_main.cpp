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
  bool verboseAll = verbose;
  verbose = (verbose && (myImageID == 0));
  
  // start the testing
	if(verbose) outputStartMessage("Directory");
  int ierr = 0;
  
  // call the actual test routines
  ierr += unitTests<int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<unsigned int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Directory", (ierr == 0));
  if(ierr != 0) cout << "[Image " << myImageID << "] Return Value = " << ierr << endl;
	return(ierr);
}

//======================================================================
template <typename OrdinalType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages) {
  std::string className = "Directory<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

  OrdinalType const nME = intToOrdinal<OrdinalType>(5);
  OrdinalType const negOne = Teuchos::OrdinalTraits<OrdinalType>::zero() - Teuchos::OrdinalTraits<OrdinalType>::one();
  
  // create platform and elementspace needed for directory construction
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
#endif
  
  // fixtures
  Tpetra::ElementSpace<OrdinalType> elementspace(nME, negOne, Teuchos::OrdinalTraits<OrdinalType>::zero(), platform);
  std::vector<OrdinalType> GIDs = elementspace.getMyGlobalElements();
  std::vector<OrdinalType> imageIDs;
  std::vector<OrdinalType> LIDs;

  // constructor
  if(verbose) cout << "Calling constructor..." << endl;
  Tpetra::BasicDirectory<OrdinalType> directory(elementspace);

  // copy constructor
  if(verbose) cout << "Calling copy constructor..." << endl;
  Tpetra::BasicDirectory<OrdinalType> directory2(directory);

  // getDirectoryEntries - imageIDs only
  if(verbose) cout << "Calling getDirectoryEntries(#1)..." << endl;
  directory.getDirectoryEntries(GIDs, imageIDs);

  // getDirectoryEntries - imageIDs and LIDs
  if(verbose) cout << "Calling getDirectoryEntries(#2)..." << endl;
  directory.getDirectoryEntries(GIDs, imageIDs, LIDs);
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  // no tests yet
  
	// ======================================================================
	// finish up
	// ======================================================================
  
  comm.barrier();
	if(verbose) {
		if(returnierr == 0)
      outputHeading("Unit tests for " + className + " passed.");
		else
      outputHeading("Unit tests for " + className + " failed.");
  }
	return(returnierr);
}
