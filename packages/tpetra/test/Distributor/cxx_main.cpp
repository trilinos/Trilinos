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
#include "Tpetra_Distributor.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template<typename PacketType, typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

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
  MPI_Barrier(MPI_COMM_WORLD);
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = (verbose && (myImageID == 0));
  
  // start the testing
	if(verbose) outputStartMessage("Distributor");
  int ierr = 0;

  // call the actual test routines
	ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<double, int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Distributor", (ierr == 0));
	return(ierr);
}

//======================================================================
template<typename PacketType, typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  std::string className = "Distributor<" + Tpetra::PacketTraits<PacketType>::name() + "," + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;
  
  // ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, PacketType> platform(MPI_COMM_WORLD);
  Tpetra::MpiComm<PacketType, OrdinalType> comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, PacketType> platform;
  Tpetra::SerialComm<PacketType, OrdinalType> comm;
  ierr = returnierr + ierr; // dummy usage of ierr so gcc doesn't complain about unused variable (this is only needed in serial mode)
#endif // TPETRA_MPI

  // platform constructor
  if(verbose) cout << "Calling Platform.createDistributor()..." << endl;
  comm.barrier();
  Teuchos::RefCountPtr< Tpetra::Distributor<OrdinalType> > distributorS = platform.createDistributor(); // distributor for createFromSends
  Teuchos::RefCountPtr< Tpetra::Distributor<OrdinalType> > distributorR = platform.createDistributor(); // distributor for createFromReceives

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

#ifdef TPETRA_MPI // Only do rest of testing if not in a serial build
  
  // fixtures
  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const length = intToOrdinal<OrdinalType>(numImages);
  OrdinalType const invalid = intToOrdinal<OrdinalType>(-99);
  OrdinalType numExportIDs;
  std::vector<OrdinalType> exportImageIDs;
  OrdinalType numRemoteIDs;

  // ========================================
  // test createFromSends
  // ========================================
  if(verbose) cout << "Testing createFromSends... ";
  numExportIDs = length; // send one element to each image, including ourself
  numRemoteIDs = zero; // in MPI, we should be receiving at least one

  // fill exportImageIDs with {0, 1, 2, ... numImages-1}
  for(OrdinalType i = zero; i < length; i++) 
    exportImageIDs.push_back(i);

  distributorS->createFromSends(numExportIDs, exportImageIDs, true, numRemoteIDs);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "exportImageIDs: " + toString(exportImageIDs));
    outputData(myImageID, numImages, "numRemoteIDs: " + toString(numRemoteIDs));
    if(verbose) cout << "[  All  ] Expected: " << numExportIDs << endl; // should be same on all images
    if(verbose) cout << "createFromSends test ";
  }
  if(numRemoteIDs != numExportIDs) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
  //distributorS->printInfo(cout);

  
  // ========================================
  // test createFromRecvs
  // ========================================
  if(verbose) cout << "Testing createFromRecvs... ";
  comm.barrier();
  numRemoteIDs = numExportIDs; // same as numExportIDs above (numImages)

  // fill remoteGIDs with row from generator: (0, myImageID), to (numImages-1, myImageID)
  std::vector<OrdinalType> remoteGIDs;
  for(OrdinalType i = zero; i < length; i++) {
    remoteGIDs.push_back(generateValue(i, intToOrdinal<OrdinalType>(myImageID)));
  }

  std::vector<OrdinalType> remoteImageIDs(exportImageIDs); // same as exportImageIDs array above {0,1,2,...(numImages-1)}
  numExportIDs = zero; // in MPI, we should be sending at least one
  std::vector<OrdinalType> exportGIDs(length, invalid);
  exportImageIDs = exportGIDs; // fill with same thing

  distributorR->createFromRecvs(numRemoteIDs, remoteGIDs, remoteImageIDs, true, numExportIDs, exportGIDs, exportImageIDs);
  std::vector<OrdinalType> expectedGIDs;
  generateColumn(expectedGIDs, myImageID, numImages);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "remoteGIDs: " + toString(remoteGIDs));
    outputData(myImageID, numImages, "numExportIDs: " + toString(numExportIDs) + " Expected: " + toString(numRemoteIDs));
    outputData(myImageID, numImages, "exportGIDs: " + toString(exportGIDs) + " Expected: " + toString(expectedGIDs));
    outputData(myImageID, numImages, "exportImageIDs: " + toString(exportImageIDs) + " Expected: " + toString(remoteImageIDs));
    if(verbose) cout << "createFromRecvs test ";
  }
  if(numExportIDs != numRemoteIDs)
    ierr++;
  if(exportGIDs != expectedGIDs)
    ierr++;
  if(exportImageIDs != remoteImageIDs)
    ierr++;
  if(ierr != 0) {
    ierr = 1;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;


  // ========================================
  // test doPosts
  // ========================================
  if(verbose) cout << "Testing doPostsAndWaits... ";
  comm.barrier();

  char* exportObjs = 0;
  OrdinalType const objectSize = Tpetra::PacketTraits<PacketType>::packetSize();
  OrdinalType lenImportObjs = invalid;
  char* importObjs = 0;
  std::vector<PacketType> exports;
  generateColumn(exports, myImageID, numImages);
  exportObjs = reinterpret_cast<char*>(&exports.front());
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "exports: " + toString(exports));
    outputData(myImageID, numImages, "objectSize: " + toString(objectSize));
  }
  distributorS->doPostsAndWaits(exportObjs, objectSize, lenImportObjs, importObjs);
  if(debug) {
    outputData(myImageID, numImages, "lenImportObjs: " + toString(lenImportObjs));
    PacketType* imports = reinterpret_cast<PacketType*>(importObjs);
    std::string importString = "{";
    for(OrdinalType i = 0; i < (lenImportObjs / Tpetra::PacketTraits<PacketType>::packetSize()); i++) {
      importString += toString(*imports++);
      importString += " ";
    }
    importString += "}";
    outputData(myImageID, numImages, "importObjs: " + importString);
    if(verbose) cout << "doPostsAndWaits test ";
  }

  if(ierr != 0) {
    ierr = 1;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

#endif // TPETRA_MPI

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
