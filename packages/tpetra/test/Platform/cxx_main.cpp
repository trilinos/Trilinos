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
#include <Teuchos_RefCountPtr.hpp>
#include "Tpetra_OmniPlatform.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template <typename PacketType, typename OrdinalType>
int omniTest(bool verbose, bool debug, int myImageID, int numImages);
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);
template <typename OrdinalType, typename ScalarType>
void codeCoverage(bool verbose, bool debug, int myImageID, int numImages);

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
	if(verbose) outputStartMessage("Platform");
	int ierr = 0;

	//omniTest<int, int>(verbose, debug, myImageID, numImages);
  
	// call the actual test routines
	ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);
	ierr += unitTests<int, complex<double> >(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
	MPI_Finalize();
#endif
	if(verbose) outputEndMessage("Platform", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename PacketType, typename OrdinalType>
int omniTest(bool verbose, bool debug, int myImageID, int numImages) {

#ifndef TPETRA_MPI
	// create Serial OmniPlatform if MPI is not enabled
	if(verbose) 
		cout << "Creating Serial OmniPlatform..." << endl;
	Tpetra::OmniPlatform op1;
	if(debug) 
		cout << op1 << endl;

	// assert that size and rank are correct
	int op1_mIID = op1.getMyImageID();
	int op1_nI = op1.getNumImages();
	if(op1_mIID != myImageID)
		cout << "** myImageID = " << myImageID << ", op1_mIID = " << op1_mIID << endl;
	if(op1_nI != numImages)
	cout << "** numImages = " << numImages << ", op1_nI = " << op1_nI << endl;

	// create the Comms
	if(verbose) 
		cout << "Creating op1's Comm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<PacketType, OrdinalType> > op1_comm;
	op1_comm = op1.createComm<PacketType, OrdinalType>();
	if(debug)
		op1_comm->printInfo(cout);
	if(verbose) 
		cout << "Creating op1's SerialComm..." << endl;
	op1_comm = op1.createSerialComm<PacketType, OrdinalType>();
	if(debug)
		op1_comm->printInfo(cout);

#else
//#ifdef TPETRA_MPI
	// create MPI OmniPlatform if MPI is enabled
	if(verbose) 
		cout << "Creating MPI OmniPlatform..." << endl;
	Tpetra::OmniPlatform op2(MPI_COMM_WORLD);
	if(debug) {
		cout << op2 << endl;
		cout.flush();
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// assert that size and rank are correct
	int op2_mIID = op2.getMyImageID();
	int op2_nI = op2.getNumImages();
	if(op2_mIID != myImageID)
		cout << "** myImageID = " << myImageID << ", op2_mIID = " << op2_mIID << endl;
	if(op2_nI != numImages)
		cout << "** numImages = " << numImages << ", op2_nI = " << op2_nI << endl;

	// create the Comms
	if(verbose) 
		cout << "Creating op2's Comm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<PacketType, OrdinalType> > op2_comm;
	op2_comm = op2.createComm<PacketType, OrdinalType>();
	if(debug)
		op2_comm->printInfo(cout);
	if(verbose) 
		cout << "Creating op2's SerialComm..." << endl;
	op2_comm = op2.createSerialComm<PacketType, OrdinalType>();
	if(debug)
		op2_comm->printInfo(cout);
#endif

	if(verbose) cout << "Finished OmniPlatform testing." << endl;

	return(0);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	std::string className = "Platform<" + Teuchos::OrdinalTraits<OrdinalType>::name() + "," + Teuchos::ScalarTraits<ScalarType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	int ierr = 0;
	int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
	codeCoverage<OrdinalType, ScalarType>(verbose, debug, myImageID, numImages);
	
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	if(verbose && debug) outputSubHeading("Starting actual testing section...");

#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform(MPI_COMM_WORLD);
	Tpetra::MpiComm<ScalarType, OrdinalType> comm(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform;
	Tpetra::SerialComm<ScalarType, OrdinalType> comm;
#endif

	// test getMyImageID
	if(verbose) cout << "Testing getMyImageID... ";
	int platform_myImageID = platform.getMyImageID();
	if(debug) {
		if(verbose) cout << endl;
		comm.barrier();
		cout << "[Image " << myImageID << "] getMyImageID: " << platform_myImageID << endl;
		cout << "[Image " << myImageID << "] Expected: " << myImageID << endl;
		comm.barrier();
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
	int platform_numImages = platform.getNumImages();
	if(debug) {
		if(verbose) cout << endl;
		comm.barrier();
		cout << "[Image " << myImageID << "] getNumImages: " << platform_numImages << endl;
		cout << "[Image " << myImageID << "] Expected: " << numImages << endl;
		comm.barrier();
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
template <typename OrdinalType, typename ScalarType>
void codeCoverage(bool verbose, bool debug, int myImageID, int numImages) { 
	if(verbose) outputSubHeading("Starting code coverage section...");

#ifdef TPETRA_MPI
	// default constructor
	if(verbose) cout << "MpiPlatform default constructor..." << endl;
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform(MPI_COMM_WORLD);
	// copy constructor
	if(verbose) cout << "MpiPlatform copy constructor..." << endl;
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform2(platform);
#else
	// default constructor
	if(verbose) cout << "SerialPlatform default constructor..." << endl;
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform;
	// copy constructor
	if(verbose) cout << "SerialPlatform copy constructor..." << endl;
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform2(platform);
#endif

	// clone
	if(verbose) cout << "clone..." << endl;
	Teuchos::RefCountPtr< Tpetra::Platform<OrdinalType, ScalarType> > platform3 = platform.clone();

	// createScalarComm
	if(verbose) cout << "createScalarComm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<ScalarType, OrdinalType> > comm1 = platform.createScalarComm();
  
	// createOrdinalComm
	if(verbose) cout << "createOrdinalComm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<OrdinalType, OrdinalType> > comm2 = platform.createOrdinalComm();

	// createDistributor
	if(verbose) cout << "createDistributor..." << endl;
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platformO(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformO;
#endif
	Teuchos::RefCountPtr< Tpetra::Distributor<OrdinalType> > distributor = platform.createDistributor();

	// create Directory
	if(verbose) cout << "createDirectory..." << endl;
	Tpetra::ElementSpace<OrdinalType> elementspace(10, 0, platformO);
	Teuchos::RefCountPtr< Tpetra::Directory<OrdinalType> > directory = platform.createDirectory(elementspace);
}
