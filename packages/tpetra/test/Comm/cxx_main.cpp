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
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template <typename PacketType, typename OrdinalType>
int simpleTest(bool const verbose, bool const debug, int const myImageID, int const numImages);
template <typename PacketType, typename OrdinalType>
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
	//simpleTest<int, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<double, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<complex<double>, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<complex<float>, int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Comm", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename PacketType, typename OrdinalType>
int simpleTest(bool const verbose, bool const debug, int const myImageID, int const numImages) {
	// compute which image to send to/from
	int destImageID = myImageID - 1;
	int sourceImageID = myImageID + 1;
	bool sending = (myImageID % 2 != 0);
	bool receiving = (!sending && sourceImageID < numImages);

	OrdinalType const numValues = intToOrdinal<OrdinalType>(4); // magic number
	int const tag = 26001; // another magic number
#ifdef TPETRA_MPI
	if(sending) {
		// initialize values we'll send
		std::vector<PacketType> myVals(numValues);
		generateColumn(myVals, myImageID, numValues);

		if(debug) 
			cout << "[Image " << myImageID << "] sending to IID " << destImageID << ", values = "
				 << Tpetra::toString(myVals) << endl;

		int ierr = MPI_Send(&myVals[0], numValues, Tpetra::MpiTraits<PacketType>::datatype(), 
							destImageID, tag, MPI_COMM_WORLD);
		if(ierr != 0)
			cerr << "[Image " << myImageID << "] Error in MPI_Send, returned " << ierr << "." << endl;
	}

	if(receiving) {
		// initialize expected values, and buffer
		// for receiving values
		std::vector<PacketType> received(numValues);
		std::vector<PacketType> expected(numValues);
		generateColumn(expected, sourceImageID, numValues);

		MPI_Status status;
		int ierr = MPI_Recv(&received[0], numValues, Tpetra::MpiTraits<PacketType>::datatype(), 
							sourceImageID, tag, MPI_COMM_WORLD, &status);
		if(ierr != 0)
			cerr << "[Image " << myImageID << "] Error in MPI_Recv, returned " << ierr << "." << endl;
		if(received != expected || debug)
			cout << "[Image " << myImageID << "] receiving from IID " << sourceImageID << ", expected = "
				 << Tpetra::toString(expected) << ", received = " << Tpetra::toString(received) << endl;
	}

	if(!receiving && !sending) {
		if(debug)
			cout << "[Image " << myImageID << "] Not participating." << endl;
	}

	// now try it through MpiComm

	Tpetra::MpiComm<PacketType, OrdinalType> comm(MPI_COMM_WORLD);
	std::vector<PacketType> sendVals(numValues);
	std::vector<PacketType> recvVals(numValues);
	
	if(myImageID == 0) { // send
		cout << "Image 0, about to send" << endl;
		comm.send(&sendVals[0], numValues, 1);
		cout << "Image 0, done" << endl;
	}
	else if(myImageID == 1) { // receive
		cout << "Image 1, about to receive" << endl;
		comm.receive(&recvVals[0], numValues, 0);
		cout << "Image 1, done" << endl;
	}
	else { // do nothing
		cout << "Image " << myImageID << ", doing nothing" << endl;
	}
	

#endif
	return(0);
}

//======================================================================
template <typename PacketType, typename OrdinalType>
int unitTests(bool const verbose, bool const debug, int const myImageID, int const numImages) {
	std::string className = "Comm<" + Tpetra::PacketTraits<PacketType>::name() + "," + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	int ierr = 0;
	int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  
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
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================
  
	// fixtures
	int const root = 0; // root image for broadcast
	OrdinalType const length = 5; // length of arrays used for Comm operations
	std::vector<PacketType> myVals(length);
	std::vector<PacketType> allVals(length);
	std::vector<PacketType> expected(length);

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
