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
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_Version.hpp"
#include "Tpetra_ElementSpace.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int rank, int size);

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
	if(verbose) {
    cout << "****************************************\n" 
    << "Starting PlatformTest..." << endl
    << Tpetra::Tpetra_Version() << endl
    << "****************************************\n";
  }
  int ierr = 0;
  
  // call the actual test routines
	ierr += unitTests<int, double>(verbose, debug, rank, size);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) {
    cout << "**************************************** " << endl;
		if(ierr == 0)
			cout << "Platform test passed." << endl;
		else
			cout << "Platform test failed." << endl;
    cout << "**************************************** " << endl;
  }
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int rank, int size) {
  if(verbose) cout << "Stating unit tests for Platform<" << Teuchos::OrdinalTraits<OrdinalType>::name() 
                   << ", " << Teuchos::ScalarTraits<ScalarType>::name() << "> " << endl;
  int ierr = 0;
  int returnierr = 0;
#ifdef TPETRA_MPI
  Tpetra::MpiComm<ScalarType, OrdinalType> comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<ScalarType, OrdinalType> comm;
#endif

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  if(verbose) cout << "Starting code coverage section..." << endl;
#ifdef TPETRA_MPI
  // default constructor
	if(verbose) cout << "MpiPlatform default constructor..." << endl;
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform(MPI_COMM_WORLD);
  if(debug) {comm.barrier(); cout << platform; comm.barrier();}
  // copy constructor
  if(verbose) cout << "MpiPlatform copy constructor..." << endl;
  Tpetra::MpiPlatform<OrdinalType, ScalarType> platform2(platform);
  if(debug) {comm.barrier(); cout << platform2; comm.barrier();}
#else
  // default constructor
	if(verbose) cout << "SerialPlatform default constructor..." << endl;
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform;
  if(debug) {comm.barrier(); cout << platform; comm.barrier();}
  // copy constructor
  if(verbose) cout << "SerialPlatform copy constructor..." << endl;
  Tpetra::SerialPlatform<OrdinalType, ScalarType> platform2(platform);
  if(debug) {comm.barrier(); cout << platform2; comm.barrier();}
#endif
  // only the constructors needed the conditionals
  // the rest of the testing is MPI/Serial independent

  // clone
  if(verbose) cout << "clone..." << endl;
  Teuchos::RefCountPtr< Tpetra::Platform<OrdinalType, ScalarType> > platform3 = platform.clone();
  if(debug) {comm.barrier(); platform3->printInfo(cout); comm.barrier();}

  // create OrdinalComm and ScalarComm
	if(verbose) cout << "createScalarComm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<ScalarType, OrdinalType> > comm1 = platform.createScalarComm();
  if(debug) {comm.barrier(); comm1->printInfo(cout); comm.barrier();}
  if(verbose) cout << "createOrdinalComm..." << endl;
	Teuchos::RefCountPtr< Tpetra::Comm<OrdinalType, OrdinalType> > comm2 = platform.createOrdinalComm();
	if(debug) {comm.barrier(); comm2->printInfo(cout); comm.barrier();}

  // create OrdinalDistributor and ScalarDistributor
	if(verbose) cout << "createScalarDistributor..." << endl;
	Teuchos::RefCountPtr< Tpetra::Distributor<ScalarType, OrdinalType> > distributor1 = platform.createScalarDistributor();
  if(debug) {comm.barrier(); distributor1->printInfo(cout); comm.barrier();}
  if(verbose) cout << "createOrdinalDistributor..." << endl;
	Teuchos::RefCountPtr< Tpetra::Distributor<OrdinalType, OrdinalType> > distributor2 = platform.createOrdinalDistributor();
	if(debug) {comm.barrier(); distributor2->printInfo(cout); comm.barrier();}

  // create Directory
  if(verbose) cout << "createDirectory..." << endl;
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD); // I lied about the rest being independent
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;              // But everything else is, I promise.
#endif
  Tpetra::ElementSpace<OrdinalType> elementspace(10, 0, platformE);
  Teuchos::RefCountPtr< Tpetra::Directory<OrdinalType> > directory = platform.createDirectory(elementspace);
  // directory doesn't have a printInfo method
	
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  if(verbose) cout << "Starting actual testing section..." << endl;

  // test getMyImageID
  if(verbose) cout << "Testing getMyImageID... ";
  int platform_rank = platform.getMyImageID();
  if(platform_rank != rank) {
    if(verbose) cout << "Failed" << endl;
    if(debug) cout << "getMyImageID returned " << platform_rank << ", should be " << rank << "." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "Passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // test getNumImages
  if(verbose) cout << "Testing getNumImages... ";
  int platform_size = platform.getNumImages();
  if(platform_size != size) {
    if(verbose) cout << "Failed" << endl;
    if(debug) cout << "getNumImages returned " << platform_size << ", should be " << size << "." << endl;
    ierr++;
  }
  else
    if(verbose) cout << "Passed" << endl;
  returnierr += ierr;
  ierr = 0;

	// ======================================================================
	// finish up
	// ======================================================================
  
	if(verbose) {
    cout << "PlatformTest <" << Teuchos::OrdinalTraits<OrdinalType>::name() 
         << ", " << Teuchos::ScalarTraits<ScalarType>::name() << "> ";
		if(returnierr == 0)
			cout << "passed." << endl;
		else
			cout << "failed." << endl;
  }
	return(returnierr);
}
