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
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#include "Tpetra_MpiTraits.hpp"
#endif // TPETRA_MPI
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"

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
	if(verbose) outputStartMessage("Comm");
	int ierr = 0;

	//mpiBreakpoint(myImageID);

	// call the actual test routines
	ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, complex<double> >(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, complex<float> >(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Comm", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages) {
	std::string className = "Comm<" + Teuchos::OrdinalTraits<OrdinalType>::name() + "," + Teuchos::ScalarTraits<ScalarType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	int ierr = 0;
	int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  
#ifdef TPETRA_MPI
	// default constructor
	if(verbose) cout << "Calling MpiComm default constructor..." << endl;
	Tpetra::MpiComm<OrdinalType, ScalarType> comm(MPI_COMM_WORLD);
	// copy constructor
	if(verbose) cout << "Calling MpiComm copy constructor..." << endl;
	Tpetra::MpiComm<OrdinalType, ScalarType> comm2(comm);
#else
	// default constructor
	if(verbose) cout << "Calling SerialComm default constructor..." << endl;
	Tpetra::SerialComm<OrdinalType, ScalarType> comm;
	// copy constructor
	if(verbose) cout << "Calling SerialComm copy constructor..." << endl;
	Tpetra::SerialComm<OrdinalType, ScalarType> comm2(comm);
#endif
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================
  
	// fixtures
	int const root = 0; // root image for broadcast
	OrdinalType const length = 5; // length of arrays used for Comm operations
	std::vector<ScalarType> myVals(length);
	std::vector<ScalarType> allVals(length);
	std::vector<ScalarType> expected(length);

	// test getMyImageID
	if(verbose) cout << "Testing getMyImageID... ";
	int platform_myImageID = comm.getMyImageID();
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "getMyImageID: " + Tpetra::toString(platform_myImageID));
		if(verbose) cout << "getMyImageID test ";
	}
	if(platform_myImageID != myImageID) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;

	// test getNumImages
	if(verbose) cout << "Testing getNumImages... ";
	int platform_numImages = comm.getNumImages();
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "getNumImages: " + Tpetra::toString(platform_numImages));
		if(verbose) cout << "getNumImages test ";
	}
	if(platform_numImages != numImages) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// test broadcast
	if(verbose) cout << "Testing broadcast... ";
	generateColumn(myVals, myImageID, length); // set myVals
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "Values prior to broadcast: " + Tpetra::toString(myVals));
	}
	comm.broadcast(&myVals.front(), length, root);
	generateColumn(expected, root, length);
	if(debug) {
		outputData(myImageID, numImages, "Values after broadcast:    " + Tpetra::toString(myVals));
		if(verbose) cout << "[  All  ] Expected values:           " << Tpetra::toString(expected) << endl;
		if(verbose) cout << "Broadcast test ";
	}
	if(myVals != expected) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
  
	// test gatherAll
	if(verbose) cout << "Testing gatherAll... ";
	generateColumn(myVals, myImageID, length); // reset myVals
	allVals.resize(length * numImages); // allVals needs to be larger for this test
	comm.gatherAll(&myVals.front(), &allVals.front(), length);
	generateMultipleColumns(expected, 0, (numImages-1), length);
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "myVals:   " + Tpetra::toString(myVals));
		outputData(myImageID, numImages, "allVals:  " + Tpetra::toString(allVals));
		if(verbose) cout << "[  All  ] Expected: " << Tpetra::toString(expected) << endl;
		if(verbose) cout << "GatherAll test ";
	}
	if(allVals != expected) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// test sumAll
	if(verbose) cout << "Testing sumAll... ";
	allVals.resize(length); // resize allVals back down
	comm.sumAll(&myVals.front(), &allVals.front(), length);
	generateRowSums(expected, 0, numImages-1, length);
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "localSums:  " + Tpetra::toString(myVals));
		outputData(myImageID, numImages, "globalSums: " + Tpetra::toString(allVals));
		if(verbose) cout << "[  All  ] Expected:   " << Tpetra::toString(expected) << endl;
		if(verbose) cout << "SumAll test ";
	}
	if(allVals != expected) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;

	// test sumAllAndScatter
	if(verbose) cout << "Testing sumAllAndScatter... ";
	generateColumn(myVals, myImageID, numImages);
	std::vector<int> recvCounts(numImages, 1);
	ScalarType receivedSum = Teuchos::ScalarTraits<OrdinalType>::zero();
	comm.sumAllAndScatter(&myVals.front(), &receivedSum, &recvCounts.front());
	generateRowSums(expected, 0, numImages-1, numImages);
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "localSums:  " + Tpetra::toString(myVals));
		outputData(myImageID, numImages, "globalSum: " + Tpetra::toString(receivedSum));
		outputData(myImageID, numImages, "Expected: " + Tpetra::toString(receivedSum));
		if(verbose) cout << "SumAllAndScatter test ";
	}
	if(receivedSum != expected[myImageID]) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	generateColumn(myVals, myImageID, length); // reset myVals
	allVals.resize(length); // reset size of allVals

  
	// test maxAll
	if(verbose) cout << "Testing maxAll... ";
	
	comm.maxAll(&myVals.front(), &allVals.front(), length);
	generateRowMaxs(expected, 0, numImages-1, length);
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "localMaxs:  " + Tpetra::toString(myVals));
		outputData(myImageID, numImages, "globalMaxs: " + Tpetra::toString(allVals));
		if(verbose) cout << "[  All  ] Expected:   " << Tpetra::toString(expected) << endl;
		if(verbose) cout << "MaxAll test ";
	}
	if(allVals != expected) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;

	// test minAll
	if(verbose) cout << "Testing minAll... ";
	comm.minAll(&myVals.front(), &allVals.front(), length);
	generateRowMins(expected, 0, numImages-1, length);
	if(debug) {
		if(verbose) cout << endl;
		comm.barrier();
		outputData(myImageID, numImages, "localMins:  " + Tpetra::toString(myVals));
		outputData(myImageID, numImages, "globalMins: " + Tpetra::toString(allVals));
		if(verbose) cout << "[  All  ] Expected:   " << Tpetra::toString(expected) << endl;
		comm.barrier();
		if(verbose) cout << "MinAll test ";
	}
	if(allVals != expected) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
  
	// test scanSum
	if(verbose) cout << "Testing scanSum... ";
	comm.scanSum(&myVals.front(), &allVals.front(), length);
	generateRowSums(expected, 0, myImageID, length);
	if(debug) {
		if(verbose) cout << endl;
		outputData(myImageID, numImages, "localScanSums:  " + Tpetra::toString(myVals));
		outputData(myImageID, numImages, "globalScanSums: " + Tpetra::toString(allVals));
		outputData(myImageID, numImages, "Expected:  " + Tpetra::toString(expected));
		if(verbose) cout << "ScanSum test ";
	}
	if(allVals != expected) {
		if(verbose) cout << "failed" << endl;
		ierr++;
	}
	else
		if(verbose) cout << "passed" << endl;
	returnierr += ierr;
	ierr = 0;
	

	#ifdef TPETRA_MPI

	// test send and receive
	cout.flush();
	comm.barrier();

	if(verbose) cout << "Testing send and receive...";
	if(debug && verbose) cout << endl;
	// pair up images - image on right sends, image on left receives
	// if there's an odd number of images, the last one doesn't participate
	// we'll send the myVals array we've been using in the previous tests
	int destImageID = myImageID - 1;
	int sourceImageID = myImageID + 1;
	bool sending = (myImageID % 2 != 0);
	bool receiving = (!sending && sourceImageID < numImages);

	if(sending) {
		generateColumn(myVals, myImageID, length);
		if(debug) {
			cout << "[Image " << myImageID << "] sending to IID " << destImageID << ", myVals = "
				 << Tpetra::toString(myVals) << endl;
		}
		comm.send(&myVals.front(), length, destImageID);
	}
	else if(receiving) {
		generateColumn(expected, sourceImageID, length);
		if(debug) {
			cout << "[Image " << myImageID << "] receiving from IID " << sourceImageID << ", expected = "
				 << Tpetra::toString(expected) << endl;
		}

		comm.receive(&allVals.front(), length, sourceImageID);

		if(debug)
			cout << "[Image " << myImageID << "] received = " << Tpetra::toString(allVals) << endl;
		if(allVals != expected)
			ierr++;

	}
	else {
		if(debug) {
			cout << "[Image " << myImageID << "] not participating" << endl;
		}
	}

	cout.flush();
	comm.barrier();

	// Note that only receivers do any checking. Senders & nonparticipants will always pass.
	if(verbose) {
		if(debug) 
			cout << "Send and Receive test ";
		if(ierr == 0)
			cout << "passed" << endl;
		else
			cout << "failed" << endl;
	}
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
