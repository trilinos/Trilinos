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

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>
#include "Tpetra_PacketTraits.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Util.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

void outputHeader(ostream& os, std::string message);

template <typename PacketType, typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size);

template <typename T>
T generateValue(int const x, int const y);
template <typename T>
void generateColumn(Teuchos::Array<T>& vector, int const x, int const length);
template <typename T>
void generateMultipleColumns(Teuchos::Array<T>& vector, int const firstx, int const lastx, int const length);
void testGenerator();

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
  if(verbose) cout << "MPI Startup: Image " << rank << " of " << size << " is alive." << endl;
  MPI_Barrier(MPI_COMM_WORLD);
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0, and verboseAll to have the original verbose setting
  bool verboseAll = verbose;
  verbose = (verbose && (rank == 0));
  
  // start the testing
	if(verbose) outputHeader(cout, "Starting CommTest...\n" + Tpetra::Tpetra_Version());
  int ierr = 0;
  
  // call the actual test routines
	//ierr += unitTests<int, int>(verbose, debug, rank, size);
	ierr += unitTests<double, int>(verbose, debug, rank, size);

  //testGenerator<int>();
  //testGenerator<float>();
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) {
		if(ierr == 0)
			outputHeader(cout, "Comm test passed.");
		else
			outputHeader(cout, "Comm test failed.");
  }
	return(ierr);
}

//======================================================================
template <typename PacketType, typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size) {
  if(verbose) cout << "Stating unit tests for Comm<" << Tpetra::PacketTraits<PacketType>::name() 
                   << ", " << Teuchos::OrdinalTraits<OrdinalType>::name() << "> " << endl;
  int ierr = 0;
  int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  if(verbose) cout << "Starting code coverage section..." << endl;
#ifdef TPETRA_MPI
  // default constructor
	if(verbose) cout << "MpiComm default constructor..." << endl;
	Tpetra::MpiComm<PacketType, OrdinalType> comm(MPI_COMM_WORLD);
  if(debug) {comm.barrier(); cout << comm; comm.barrier();}
  // copy constructor
  if(verbose) cout << "MpiComm copy constructor..." << endl;
  Tpetra::MpiComm<PacketType, OrdinalType> comm2(comm);
  if(debug) {comm.barrier(); cout << comm2; comm.barrier();}
#else
  // default constructor
	if(verbose) cout << "SerialComm default constructor..." << endl;
	Tpetra::SerialComm<PacketType, OrdinalType> comm;
  if(debug) {comm.barrier(); cout << comm; comm.barrier();}
  // copy constructor
  if(verbose) cout << "SerialComm copy constructor..." << endl;
  Tpetra::SerialComm<PacketType, OrdinalType> comm2(comm);
  if(debug) {comm.barrier(); cout << comm2; comm.barrier();}
#endif
  // only the constructors needed the conditionals
  // the rest of the testing is MPI/Serial independent

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  if(verbose) cout << "Starting actual testing section... " << endl;

  OrdinalType const length = 5;
  int const root = 0; // root image for broadcast

  // give each image a column from the generator
  // do not modify these values - they are used by all tests
  Teuchos::Array<PacketType> myVals;
  generateColumn(myVals, rank, length);

  // test broadcast
  if(verbose) cout << "Testing broadcast..." << endl;
  Teuchos::Array<PacketType> rootVals;
  generateColumn(rootVals, root, length);
  if(debug) {
    comm.barrier();
    cout << "[Image " << rank << "] Values prior to broadcast: " << myVals.toString() << endl;
  }
  comm.broadcast(&myVals[0], length, root);
  if(debug) {
    cout << "[Image " << rank << "] Values after broadcast:    " << myVals.toString() << endl;
    cout << "[Image " << rank << "] Expected values:           " << rootVals.toString() << endl;
    comm.barrier();
  }
  if(myVals != rootVals) {
    if(verbose) cout << "Broadcast test failed." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "Broadcast test passed." << endl;
  returnierr += ierr;
  ierr = 0;

  // test gatherAll
  if(verbose) cout << "Testing gatherAll..." << endl;
  generateColumn(myVals, rank, length);
  Teuchos::Array<PacketType> allVals(size * length);
  comm.gatherAll(&myVals[0], &allVals[0], length);
  Teuchos::Array<PacketType> expectedAllVals;
  generateMultipleColumns(expectedAllVals, 0, (size-1), length);
  if(debug) {
    comm.barrier();
    cout << "[Image " << rank << "] myVals:   " << myVals.toString() << endl;
    cout << "[Image " << rank << "] allVals:  " << allVals.toString() << endl;
    cout << "[Image " << rank << "] Expected: " << expectedAllVals.toString() << endl;
    comm.barrier();
  }
  if(allVals != expectedAllVals) {
    if(verbose) cout << "GatherAll test failed." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "GatherAll test passed." << endl;
  returnierr += ierr;
  ierr = 0;

  // test sumAll
  if(verbose) cout << "Testing sumAll..." << endl;
  int indexToUse = 3;
  PacketType localSum = myVals[indexToUse];
  PacketType globalSum = Teuchos::ScalarTraits<PacketType>::zero();
  comm.sumAll(&localSum, &globalSum, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedGlobalSum = Teuchos::ScalarTraits<PacketType>::zero();
  for(int x = 0; x < size; x++)
    expectedGlobalSum += generateValue<PacketType>(x, indexToUse);
  if(debug) {
    comm.barrier();
    cout << "[Image " << rank << "] localSum:  " << localSum << endl;
    cout << "[Image " << rank << "] globalSum: " << globalSum << endl;
    cout << "[Image " << rank << "] Expected:  " << expectedGlobalSum << endl;
    comm.barrier();
  }
  if(globalSum != expectedGlobalSum) {
    if(verbose) cout << "SumAll test failed." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "SumAll test passed." << endl;
  returnierr += ierr;
  ierr = 0;

  // test maxAll
  if(verbose) cout << "Testing maxAll..." << endl;
  PacketType localMax = myVals[length-1]; // dependent on generator ordering
  PacketType globalMax = Teuchos::ScalarTraits<PacketType>::zero();
  comm.maxAll(&localMax, &globalMax, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedGlobalMax = generateValue<PacketType>(size-1, length-1); // dependent on generator ordering
  if(debug) {
    comm.barrier();
    cout << "[Image " << rank << "] localMax:  " << localMax << endl;
    cout << "[Image " << rank << "] globalMax: " << globalMax << endl;
    cout << "[Image " << rank << "] Expected:  " << expectedGlobalMax << endl;
    comm.barrier();
  }
  if(globalMax != expectedGlobalMax) {
    if(verbose) cout << "MaxAll test failed." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "MaxAll test passed." << endl;
  returnierr += ierr;
  ierr = 0;

  // test minAll
  if(verbose) cout << "Testing minAll..." << endl;
  indexToUse = 1;
  PacketType localMin = myVals[indexToUse]; // might not be true min
  PacketType globalMin = Teuchos::ScalarTraits<PacketType>::zero();
  comm.minAll(&localMin, &globalMin, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedGlobalMin = generateValue<PacketType>(0, indexToUse); // dependent on generator ordering
  if(debug) {
    comm.barrier();
    cout << "[Image " << rank << "] localMin:  " << localMin << endl;
    cout << "[Image " << rank << "] globalMin: " << globalMin << endl;
    cout << "[Image " << rank << "] Expected:  " << expectedGlobalMin << endl;
    comm.barrier();
  }
  if(globalMin != expectedGlobalMin) {
    if(verbose) cout << "MinAll test failed." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "MinAll test passed." << endl;
  returnierr += ierr;
  ierr = 0;

  // test scanSum
  if(verbose) cout << "Testing scanSum..." << endl;
  indexToUse = 2;
  PacketType localScanSum = myVals[indexToUse];
  PacketType globalScanSum = Teuchos::ScalarTraits<PacketType>::zero();
  comm.scanSum(&localScanSum, &globalScanSum, Teuchos::OrdinalTraits<OrdinalType>::one());
  PacketType expectedScanSum = Teuchos::ScalarTraits<PacketType>::zero();
  for(int x = 0; x <= rank; x++)
    expectedScanSum += generateValue<PacketType>(x, indexToUse);
  if(debug) {
    comm.barrier();
    cout << "[Image " << rank << "] localScanSum:  " << localScanSum << endl;
    cout << "[Image " << rank << "] globalScanSum: " << globalScanSum << endl;
    cout << "[Image " << rank << "] Expected:      " << expectedScanSum << endl;
    comm.barrier();
  }
  if(globalScanSum != expectedScanSum) {
    if(verbose) cout << "ScanSum test failed." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "ScanSum test passed." << endl;
  returnierr += ierr;
  ierr = 0;

	// ======================================================================
	// finish up
	// ======================================================================
  
  comm.barrier();
	if(verbose) {
    cout << "CommTest <" << Tpetra::PacketTraits<PacketType>::name() 
         << ", " << Teuchos::OrdinalTraits<OrdinalType>::name() << "> ";
		if(returnierr == 0)
			cout << "passed." << endl;
		else
			cout << "failed." << endl;
  }
	return(returnierr);
}


//======================================================================
//======================================================================
void outputHeader(ostream& os, std::string message) {
  os << "****************************************\n"
     << message << endl
     << "****************************************\n";
}

//======================================================================
template <typename T>
T generateValue(int const x, int const y) {
  T const two = Teuchos::ScalarTraits<T>::one() + Teuchos::ScalarTraits<T>::one();
  // formula for z(x,y) = 0.5(x^2 + y^2 + 3x + y) + xy
  return(((x*x + y*y + x+x+x + y) / two) + (x*y));
}

//======================================================================
template <typename T>
void generateColumn(Teuchos::Array<T>& vector, int const x, int const length) {
  vector.resize(length);
  for(int y = 0; y < length; y++)
    vector[y] = generateValue<T>(x, y);
}

//======================================================================
template <typename T>
void generateMultipleColumns(Teuchos::Array<T>& vector, int const firstx, int const lastx, int const length) {
  vector.resize(length * (lastx - firstx + 1));
  typename Teuchos::Array<T>::iterator i = vector.begin();
  for(int x = firstx; x <= lastx; x++)
    for(int y = 0; y < length; y++) {
      *i = generateValue<T>(x, y);
      i++;
    }
  /*cout << "firstx = " << firstx << ", lastx = " << lastx << ", length = " << length << endl;
  cout << "array size = " << vector.size() << endl;
  cout << "array contents = " << vector.toString() << endl;*/
}

//======================================================================
template <typename T>
void testGenerator() {
  int length(10);

  cout << "\nTesting generateValue..." << endl;
  for(int y = 0; y < length; y++) {
    for(int x = 0; x < length; x++)
      cout << generateValue<T>(x,y) << "\t";
    cout << endl;
  }

  cout << "\nTesting generateColumn..." << endl;
  Teuchos::Array< Teuchos::Array<T> > columns(length);
  for(int i = 0; i < length; i++) {
    generateColumn(columns[i], i, length);
    cout << "Column " << i << ": " << columns[i].toString() << endl;
  }
  cout << "\n";
}
