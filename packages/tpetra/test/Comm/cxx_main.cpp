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
int unitTests(bool verbose, bool debug, int rank, int size);
template <typename PacketType, typename OrdinalType>
void codeCoverage(bool verbose, bool debug, int rank, int size);

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

  int rank = 0; // assume we are on serial
  int size = 1; // if MPI, will be reset later
  
  // initialize MPI if needed
#ifdef TPETRA_MPI
  size = -1;
  rank = -1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = (verbose && (rank == 0));
  
  // start the testing
	if(verbose) outputStartMessage("Comm");
  int ierr = 0;
  
  // call the actual test routines
  ierr += unitTests<int, int>(verbose, debug, rank, size);
	ierr += unitTests<double, int>(verbose, debug, rank, size);
  ierr += unitTests<complex<double>, int>(verbose, debug, rank, size);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Comm", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename PacketType, typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size) {
  std::string className = "Comm<" + Tpetra::PacketTraits<PacketType>::name() + "," + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  codeCoverage<PacketType, OrdinalType>((verbose && debug), rank, size);

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  if(verbose && debug) outputSubHeading("Starting actual testing section...");

#ifdef TPETRA_MPI
	Tpetra::MpiComm<PacketType, OrdinalType> comm(MPI_COMM_WORLD);
#else
	Tpetra::SerialComm<PacketType, OrdinalType> comm;
#endif

  OrdinalType const length = 5;
  int const root = 0; // root image for broadcast

  // give each image a column from the generator
  // do not modify these values - they are used by all tests
  std::vector<PacketType> myVals;
  generateColumn(myVals, rank, length);

  // test broadcast
  if(verbose) cout << "Testing broadcast... ";
  std::vector<PacketType> rootVals;
  generateColumn(rootVals, root, length);
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << rank << "] Values prior to broadcast: " << myVals << endl;
  }
  comm.broadcast(&myVals.front(), length, root);
  if(debug) {
    cout << "[Image " << rank << "] Values after broadcast:    " << myVals << endl;
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
  generateColumn(myVals, rank, length);
  std::vector<PacketType> allVals(size * length);
  comm.gatherAll(&myVals.front(), &allVals.front(), length);
  std::vector<PacketType> expectedAllVals;
  generateMultipleColumns(expectedAllVals, 0, (size-1), length);
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << rank << "] myVals:   " << myVals << endl;
    cout << "[Image " << rank << "] allVals:  " << allVals << endl;
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
  for(int x = 0; x < size; x++)
    expectedGlobalSum += generateValue(intToScalar<PacketType>(x), intToScalar<PacketType>(indexToUse));
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << rank << "] localSum:  " << localSum << endl;
    cout << "[Image " << rank << "] globalSum: " << globalSum << endl;
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
  PacketType expectedGlobalMax = generateValue(intToScalar<PacketType>(size-1), intToScalar<PacketType>(length-1)); // dependent on generator ordering
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << rank << "] localMax:  " << localMax << endl;
    cout << "[Image " << rank << "] globalMax: " << globalMax << endl;
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
    cout << "[Image " << rank << "] localMin:  " << localMin << endl;
    cout << "[Image " << rank << "] globalMin: " << globalMin << endl;
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
  for(int x = 0; x <= rank; x++)
    expectedScanSum += generateValue(intToScalar<PacketType>(x), intToScalar<PacketType>(indexToUse));
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << rank << "] localScanSum:  " << localScanSum << endl;
    cout << "[Image " << rank << "] globalScanSum: " << globalScanSum << endl;
    cout << "[Image " << rank << "] Expected:      " << expectedScanSum << endl;
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

//======================================================================
template <typename PacketType, typename OrdinalType>
void codeCoverage(bool verbose, int rank, int size) { 
  if(verbose) outputSubHeading("Starting code coverage section...");

  OrdinalType length = intToOrdinal<OrdinalType>(10);
  std::vector<PacketType> myVals(length);
  std::vector<PacketType> allVals(length);
  int root = 0;

#ifdef TPETRA_MPI
  // default constructor
	if(verbose) cout << "MpiComm default constructor..." << endl;
	Tpetra::MpiComm<PacketType, OrdinalType> comm(MPI_COMM_WORLD);
  // copy constructor
  if(verbose) cout << "MpiComm copy constructor..." << endl;
  Tpetra::MpiComm<PacketType, OrdinalType> comm2(comm);
#else
  // default constructor
	if(verbose) cout << "SerialComm default constructor..." << endl;
	Tpetra::SerialComm<PacketType, OrdinalType> comm;
  // copy constructor
  if(verbose) cout << "SerialComm copy constructor..." << endl;
  Tpetra::SerialComm<PacketType, OrdinalType> comm2(comm);
#endif

  // barrier
  if(verbose) cout << "barrier..." << endl;
  comm.barrier();

  // broadcast
  if(verbose) cout << "broadcast..." << endl;
  comm.broadcast(&myVals.front(), length, root);
  
  // gatherAll
  if(verbose) cout << "gatherAll..." << endl;
  comm.gatherAll(&myVals.front(), &allVals.front(), length);
  
  // sumAll
  if(verbose) cout << "sumAll..." << endl;
  comm.sumAll(&myVals.front(), &allVals.front(), length);
  
  // maxAll
  if(verbose) cout << "maxAll..." << endl;
  comm.maxAll(&myVals.front(), &allVals.front(), length);

  // minAll
  if(verbose) cout << "minAll..." << endl;
  comm.minAll(&myVals.front(), &allVals.front(), length);

  // scanSum
  if(verbose) cout << "scanSum..." << endl;
  comm.scanSum(&myVals.front(), &allVals.front(), length);
  
  // printInfo
  //if(verbose) cout << "printInfo..." << endl;
  //ostream& blackhole = Teuchos::basic_oblackholestream<void, void>();
  //comm.printInfo(blackhole);
}
