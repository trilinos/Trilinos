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
#include "Tpetra_Import.hpp"
#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
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
	// if debug is enabled, it will still output on all nodes
	verbose = (verbose && (myImageID == 0));
  
	// start the testing
	if(verbose) outputStartMessage("ImportExport");
	int ierr = 0;

	// call the actual test routines
	ierr += unitTests<int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("ImportExport", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	std::string className = "Import<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
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
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
	Tpetra::MpiComm<OrdinalType, OrdinalType> comm(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
	Tpetra::SerialComm<OrdinalType, OrdinalType> comm;
#endif

	// create Source and Target ElementSpaces
	Tpetra::ElementSpace<OrdinalType> source(negOne, intToOrdinal<OrdinalType>(10), zero, platform);
	Tpetra::ElementSpace<OrdinalType> target(negOne, intToOrdinal<OrdinalType>(5), zero, platform);
	
	// Import constructors
	if(verbose) cout << "Import constructor..." << endl;
	Tpetra::Import<OrdinalType> importer(source, target);

	if(verbose) cout << "Import copy constructor..." << endl;
	Tpetra::Import<OrdinalType> importer2(importer);

	// attribute accessors
	if(verbose) cout << "getPermuteFromLIDs..." << endl;
	importer.getPermuteFromLIDs();

	if(verbose) cout << "getPermuteToLIDs..." << endl;
	importer.getPermuteToLIDs();

	if(verbose) cout << "getRemoteIDs..." << endl;
	importer.getRemoteLIDs();

	if(verbose) cout << "getNumExportIDs..." << endl;
	importer.getNumExportIDs();

	if(verbose) cout << "getExportLIDs..." << endl;
	importer.getExportLIDs();

	if(verbose) cout << "getExportImageIDs..." << endl;
	importer.getExportImageIDs();

	if(verbose) cout << "getSourceSpace..." << endl;
	importer.getSourceSpace();

	if(verbose) cout << "getTargetSpace..." << endl;
	importer.getTargetSpace();

	// assignment
	if(verbose) cout << "assignment operator..." << endl;
	importer2 = importer;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	// test that numSame + numPermute + numRemote = target.getNumMyElements()
	if(verbose) cout << "Testing same/permute/remote sum... ";

	OrdinalType same = importer.getNumSameIDs();
	OrdinalType permute = importer.getNumPermuteIDs();
	OrdinalType remote = importer.getNumRemoteIDs();
	OrdinalType sum = same + permute + remote;
	OrdinalType expectedSum = target.getNumMyElements();
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "NumSameIDs: " + Tpetra::toString(same));
		outputData(myImageID, numImages, "NumPermuteIDs: " + Tpetra::toString(permute));
		outputData(myImageID, numImages, "NumRemoteIDs: " + Tpetra::toString(remote));
		outputData(myImageID, numImages, "Expected Sum: " + Tpetra::toString(expectedSum));
		if(verbose) cout << "same/permute/remote sum test ";
	}
	if(sum != expectedSum) {
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
  
	comm.barrier();
	if(verbose) {
		if(returnierr == 0)
			outputHeading("Unit tests for " + className + " passed.");
		else
			outputHeading("Unit tests for " + className + " failed.");
	}
	return(returnierr);
}
