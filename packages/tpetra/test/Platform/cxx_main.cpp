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
	int mIID = -1; // dummy value
	int nI = -1; // dummy value

//#ifndef TPETRA_MPI
	// create Serial OmniPlatform if MPI is not enabled
	if(verbose) 
		cout << "Creating Serial OmniPlatform..." << endl;
	Tpetra::OmniPlatform op1;
	if(debug) 
		cout << op1 << endl;

	// create the Comms
	if(verbose) 
		cout << "Creating op1's Comm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<PacketType, OrdinalType> > op1_comm;
	op1.createComm(op1_comm);
	if(debug)
		op1_comm->printInfo(cout);
	if(verbose) 
		cout << "Creating op1's SerialComm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<PacketType, OrdinalType> > op1_serialcomm;
	op1.createSerialComm(op1_serialcomm);
	if(debug)
		op1_serialcomm->printInfo(cout);

	// assert that size and rank are correct
	mIID = op1_comm->getMyImageID();
	nI = op1_comm->getNumImages();
	if(mIID != 0)
		cout << "** op1_comm.getMyImageID() = " << mIID << ", expected 0" << endl;
	if(nI != 1)
		cout << "** op1_comm.getNumImages() = " << nI << ", expected 1" << endl;
	mIID = op1_serialcomm->getMyImageID();
	nI = op1_serialcomm->getNumImages();
	if(mIID != 0)
		cout << "** op1_serialcomm.getMyImageID() = " << mIID << ", expected 0" << endl;
	if(nI != 1)
		cout << "** op1_serialcomm.getNumImages() = " << nI << ", expected 1" << endl;

//#else
#ifdef TPETRA_MPI
	// create MPI OmniPlatform if MPI is enabled
	if(verbose) 
		cout << "Creating MPI OmniPlatform..." << endl;
	Tpetra::OmniPlatform op2(MPI_COMM_WORLD);
	if(debug) {
		cout << op2 << endl;
		cout.flush();
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// create the Comms
	if(verbose) 
		cout << "Creating op2's Comm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<PacketType, OrdinalType> > op2_comm;
	op2.createComm(op2_comm);
	if(debug)
		op2_comm->printInfo(cout);
	if(verbose) 
		cout << "Creating op2's SerialComm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<PacketType, OrdinalType> > op2_serialcomm;
	op2.createSerialComm(op2_serialcomm);
	if(debug)
		op2_serialcomm->printInfo(cout);

	// assert that size and rank are correct
	mIID = op2_comm->getMyImageID();
	nI = op2_comm->getNumImages();
	if(mIID != 0)
		cout << "** op2_comm.getMyImageID() = " << mIID << ", expected 0" << endl;
	if(nI != 1)
		cout << "** op2_comm.getNumImages() = " << nI << ", expected 1" << endl;
	mIID = op2_serialcomm->getMyImageID();
	nI = op2_serialcomm->getNumImages();
	if(mIID != 0)
		cout << "** op2_serialcomm.getMyImageID() = " << mIID << ", expected 0" << endl;
	if(nI != 1)
		cout << "** op2_serialcomm.getNumImages() = " << nI << ", expected 1" << endl;
#endif

	if(verbose) cout << "Finished OmniPlatform testing." << endl;

	return(0);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
	std::string className = "Platform<" + Teuchos::OrdinalTraits<OrdinalType>::name() + "," + Teuchos::ScalarTraits<ScalarType>::name() + ">";
	if(verbose) outputHeading("Stating unit tests for " + className);

	//int ierr = 0; // unused now that there's no actual tests
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

	// create Directory
	if(verbose) cout << "createDirectory..." << endl;
#ifdef TPETRA_MPI
	Tpetra::MpiPlatform<OrdinalType, OrdinalType> platformO(MPI_COMM_WORLD);
#else
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformO;
#endif
	Tpetra::ElementSpace<OrdinalType> elementspace(10, 0, platformO);
	Teuchos::RefCountPtr< Tpetra::Directory<OrdinalType> > directory = platform.createDirectory(elementspace);
}
