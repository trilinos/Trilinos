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
#include "Tpetra_PacketTraits.hpp"
#include "Tpetra_Util.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template <typename PacketType, typename OrdinalType>
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
	if(verbose) outputStartMessage("Comm");
  int ierr = 0;
  
  // call the actual test routines
  ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<double, int>(verbose, debug, myImageID, numImages);
  //ierr += unitTests<complex<double>, int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Comm", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename PacketType, typename OrdinalType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages) {
  std::string className = "Comm<" + Tpetra::PacketTraits<PacketType>::name() + "," + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;
  int const root = 0; // root image for broadcast
  OrdinalType const length = 5; // length of arrays used for Comm operations
  
  // give each image a column from the generator
  // do not modify these values - they are used by all tests
  std::vector<PacketType> myVals;
  generateColumn(myVals, myImageID, length);

	// ======================================================================
	// code coverage section - just call functions, no testing
	// =====================================================================
  
#ifdef TPETRA_MPI
  // default constructor
	if(verbose) cout << "Calling MpiComm default constructor..." << endl;
	Tpetra::MpiComm<PacketType, OrdinalType> comm(MPI_COMM_WORLD);
  // copy constructor
  if(verbose) cout << "Calling MpiComm copy constructor..." << endl;
  Tpetra::MpiComm<PacketType, OrdinalType> comm2(comm);
#else
  // default constructor
	if(verbose) cout << "Calling SerialComm default constructor..." << endl;
	Tpetra::SerialComm<PacketType, OrdinalType> comm;
  // copy constructor
  if(verbose) cout << "Calling SerialComm copy constructor..." << endl;
  Tpetra::SerialComm<PacketType, OrdinalType> comm2(comm);
#endif
  
  // printInfo
  //if(verbose) cout << "Calling printInfo..." << endl;
  //ostream& blackhole = Teuchos::basic_oblackholestream<void, void>();
  //comm.printInfo(blackhole);
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  // test broadcast
  if(verbose) cout << "Testing broadcast... ";
  std::vector<PacketType> rootVals;
  generateColumn(rootVals, root, length);
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << myImageID << "] Values prior to broadcast: " << myVals << endl;
  }
  comm.broadcast(&myVals.front(), length, root);
  if(debug) {
    cout << "[Image " << myImageID << "] Values after broadcast:    " << myVals << endl;
    if(verbose) cout << "[  All  ] Expected values:           " << rootVals << endl;
    comm.barrier();
    if(verbose) cout << "Broadcast test ";
  }
  if(myVals != rootVals) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
  
  // test gatherAll
  if(verbose) cout << "Testing gatherAll... ";
  generateColumn(myVals, myImageID, length);
  std::vector<PacketType> allVals(numImages * length);
  comm.gatherAll(&myVals.front(), &allVals.front(), length);
  std::vector<PacketType> expectedAllVals;
  generateMultipleColumns(expectedAllVals, 0, (numImages-1), length);
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << myImageID << "] myVals:   " << myVals << endl;
    cout << "[Image " << myImageID << "] allVals:  " << allVals << endl;
    if(verbose) cout << "[  All  ] Expected: " << expectedAllVals << endl;
    comm.barrier();
    if(verbose) cout << "GatherAll test ";
  }
  if(allVals != expectedAllVals) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
  
  // test sumAll
  if(verbose) cout << "Testing sumAll... ";
  int indexToUse = 3;
  PacketType localSum = myVals[indexToUse];
  PacketType globalSum = Teuchos::ScalarTraits<PacketType>::zero();
  comm.sumAll(&localSum, &globalSum, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedGlobalSum = Teuchos::ScalarTraits<PacketType>::zero();
  for(int x = 0; x < numImages; x++)
    expectedGlobalSum += generateValue(intToScalar<PacketType>(x), intToScalar<PacketType>(indexToUse));
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << myImageID << "] localSum:  " << localSum << endl;
    cout << "[Image " << myImageID << "] globalSum: " << globalSum << endl;
    if(verbose) cout << "[  All  ] Expected:  " << expectedGlobalSum << endl;
    comm.barrier();
    if(verbose) cout << "SumAll test ";
  }
  if(globalSum != expectedGlobalSum) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
  
  // test maxAll
  if(verbose) cout << "Testing maxAll... ";
  PacketType localMax = myVals[length-1]; // dependent on generator ordering
  PacketType globalMax = Teuchos::ScalarTraits<PacketType>::zero();
  comm.maxAll(&localMax, &globalMax, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedGlobalMax = generateValue(intToScalar<PacketType>(numImages-1), intToScalar<PacketType>(length-1)); // dependent on generator ordering
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << myImageID << "] localMax:  " << localMax << endl;
    cout << "[Image " << myImageID << "] globalMax: " << globalMax << endl;
    if(verbose) cout << "[  All  ] Expected:  " << expectedGlobalMax << endl;
    comm.barrier();
    if(verbose) cout << "MaxAll test ";
  }
  if(globalMax != expectedGlobalMax) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // test minAll
  if(verbose) cout << "Testing minAll... ";
  indexToUse = 1;
  PacketType localMin = myVals[indexToUse]; // might not be true min in whole myVals array
  PacketType globalMin = Teuchos::ScalarTraits<PacketType>::zero();
  comm.minAll(&localMin, &globalMin, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedGlobalMin = generateValue(intToScalar<PacketType>(0), intToScalar<PacketType>(indexToUse)); // dependent on generator ordering
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << myImageID << "] localMin:  " << localMin << endl;
    cout << "[Image " << myImageID << "] globalMin: " << globalMin << endl;
    if(verbose) cout << "[  All  ] Expected:  " << expectedGlobalMin << endl;
    comm.barrier();
    if(verbose) cout << "MinAll test ";
  }
  if(globalMin != expectedGlobalMin) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
  
  // test scanSum
  if(verbose) cout << "Testing scanSum... ";
  indexToUse = 2;
  PacketType localScanSum = myVals[indexToUse];
  PacketType globalScanSum = Teuchos::ScalarTraits<PacketType>::zero();
  comm.scanSum(&localScanSum, &globalScanSum, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedScanSum = Teuchos::ScalarTraits<PacketType>::zero();
  for(int x = 0; x <= myImageID; x++)
    expectedScanSum += generateValue(intToScalar<PacketType>(x), intToScalar<PacketType>(indexToUse));
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << myImageID << "] localScanSum:  " << localScanSum << endl;
    cout << "[Image " << myImageID << "] globalScanSum: " << globalScanSum << endl;
    cout << "[Image " << myImageID << "] Expected:      " << expectedScanSum << endl;
    comm.barrier();
    if(verbose) cout << "ScanSum test ";
  }
  if(globalScanSum != expectedScanSum) {
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
