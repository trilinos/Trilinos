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
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 
#ifdef HAVE_MPI
# include "Tpetra_MpiPlatform.hpp"
#else
# include "Tpetra_SerialPlatform.hpp"
#endif // HAVE_MPI
#include "Tpetra_Distributor.hpp"

using namespace Teuchos;

template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

int main(int argc, char* argv[]) {
  int myImageID = 0; // assume we are on serial
  int numImages = 1; // if MPI, will be reset later

  // initialize MPI if needed
#ifdef HAVE_MPI
  numImages = -1;
  myImageID = -1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numImages);
  MPI_Comm_rank(MPI_COMM_WORLD, &myImageID);
#endif // HAVE_MPI

  bool verbose = false;
  bool debug = false;
  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet"  ,&verbose,"Print messages and results.");
  cmdp.setOption("debug"  ,"nodebug",&debug  ,"Print debugging info.");
  cmdp.parse(argc,argv);
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = ((verbose || debug) && (myImageID == 0));


  // start the testing
  if(verbose) outputStartMessage("Distributor");
  int ierr = 0;

  // call the actual test routines
  ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
  ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);

  // finish up
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  if(verbose) outputEndMessage("Distributor", (ierr == 0));
  return(ierr);
}


//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  std::string className = "Distributor<" + OrdinalTraits<OrdinalType>::name() + "> with " + ScalarTraits<ScalarType>::name();
  if(verbose) outputHeading("Stating unit tests for " + className);

  (void)debug;
  (void)myImageID;
  (void)numImages;

  int ierr = 0;
  int returnierr = 0;

  // ======================================================================
  // code coverage section - just call functions, no testing
  // ======================================================================
#ifdef HAVE_MPI
  Tpetra::MpiPlatform<OrdinalType>    platform(rcp(new OpaqueWrapper<MPI_Comm>(MPI_COMM_WORLD)) );
#else
  Tpetra::SerialPlatform<OrdinalType> platform;
  (void)ierr;
#endif // HAVE_MPI
  RCP< Comm<OrdinalType> > comm = platform.createComm();

  // platform constructor
  if(verbose) cout << "Calling constructor..." << endl;
  comm->barrier();
  Tpetra::Distributor<OrdinalType> distributorS(comm); // distributor for createFromSends
  Tpetra::Distributor<OrdinalType> distributorR(comm); // distributor for createFromReceives

  // ======================================================================
  // actual testing section - affects return code
  // ======================================================================

#ifdef HAVE_MPI // Only do rest of testing if not in a serial build
  // fixtures
  OrdinalType const zero = OrdinalTraits<OrdinalType>::zero();
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

  distributorS.createFromSends(numExportIDs, exportImageIDs, true, numRemoteIDs);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "exportImageIDs: " + Tpetra::toString(exportImageIDs));
    outputData(myImageID, numImages, "numRemoteIDs: " + Tpetra::toString(numRemoteIDs));
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
  comm->barrier();
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

  distributorR.createFromRecvs(numRemoteIDs, remoteGIDs, remoteImageIDs, true, numExportIDs, exportGIDs, exportImageIDs);
  std::vector<OrdinalType> expectedGIDs;
  generateColumn(expectedGIDs, myImageID, numImages);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "remoteGIDs: " + Tpetra::toString(remoteGIDs));
    outputData(myImageID, numImages, "numExportIDs: " + Tpetra::toString(numExportIDs) + " Expected: " + Tpetra::toString(numRemoteIDs));
    outputData(myImageID, numImages, "exportGIDs: " + Tpetra::toString(exportGIDs) + " Expected: " + Tpetra::toString(expectedGIDs));
    outputData(myImageID, numImages, "exportImageIDs: " + Tpetra::toString(exportImageIDs) + " Expected: " + Tpetra::toString(remoteImageIDs));
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
  comm->barrier();

  OrdinalType const objectSize = OrdinalTraits<OrdinalType>::one();
  std::vector<ScalarType> imports;
  std::vector<ScalarType> exports;
  generateColumn(exports, myImageID, numImages);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "exports: " + Tpetra::toString(exports));
    outputData(myImageID, numImages, "objectSize: " + Tpetra::toString(objectSize));
  }
  // FINISH
  distributorS.doPostsAndWaits(exports, objectSize, imports);
  if(debug) {
    outputData(myImageID, numImages, "imports: " + Tpetra::toString(imports));
    if(verbose) cout << "doPostsAndWaits test: ";
  }

  if(ierr != 0) {
    ierr = 1;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
#endif // HAVE_MPI

  // ======================================================================
  // finish up
  // ======================================================================

  comm->barrier();
  if(verbose) {
    if(returnierr == 0)
      outputHeading("Unit tests for " + className + " passed.");
    else
      outputHeading("Unit tests for " + className + " failed.");
  }
  return(returnierr);
}
