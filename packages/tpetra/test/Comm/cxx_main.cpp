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

#include <iostream>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_Version.hpp"
#include "Tpetra_SerialComm.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiComm.hpp"
#endif // TPETRA_MPI

// function prototypes
template<typename PacketType, typename OrdinalType> 
int commMethods(Tpetra::Comm<PacketType, OrdinalType>& comm, bool verbose);
template<typename PacketType, typename OrdinalType> 
int checkRankAndSize(Tpetra::Comm<PacketType, OrdinalType>& comm, int rank, int size, bool verbose, bool debug);

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
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0, and verboseAll to have the original verbose setting
  bool verboseAll = verbose;
  verbose = (verbose && (rank == 0));
  
  // start the testing
	if(verbose) {
    cout << "\n****************************************\n" 
         << "Starting CommTest..." << endl
         << Tpetra::Tpetra_Version() << endl
         << "****************************************\n";
  }
  int ierr = 0;
  int returnierr = 0;
  
  // test SerialComm
  if(verbose) cout << "\n****************************************\n" 
                   << "Testing Tpetra::SerialComm\n" 
                   << "****************************************\n\n";
  if(verbose) cout << "Creating SerialComm object...";
  Tpetra::SerialComm<float, int> comm;
  if(verbose) cout << "Successful" << endl;
  Tpetra::SerialComm<float, int> serialClone(comm); // test copy constructor
  ierr += checkRankAndSize<float, int>(comm, 0, 1, verbose, debug); // for serial, rank and size are always 0 and 1
  ierr += commMethods<float, int>(comm, verbose);
  if(verbose) 
    if(ierr == 0) cout << "SerialComm testing successful." << endl;
    else cout << "SerialComm testing failed." << endl;
  returnierr += ierr;
  ierr = 0;
  
  // test MpiComm if MPI is enabled
#ifdef TPETRA_MPI
  if(verbose) cout << "\n****************************************\n" 
                   << "Testing Tpetra::MpiComm\n" 
                   << "****************************************\n\n";
  
  if(verbose) cout << "Creating MpiComm object...";
  Tpetra::MpiComm<float, int> comm2(MPI_COMM_WORLD);
  if(verbose) cout << "Successful" << endl;
  Tpetra::MpiComm<float, int> mpiClone(comm2); // test copy constructor
  ierr += checkRankAndSize<float, int>(comm2, rank, size, verbose, debug);
  ierr += commMethods<float, int>(comm2, verbose);
  if(verbose) 
    if(ierr == 0) cout << "MpiComm testing successful." << endl;
    else cout << "MpiComm testing failed." << endl;
  returnierr += ierr;
  ierr = 0;
  
  MPI_Finalize();
#endif // TPETRA_MPI
  
  // finish up
	if(verbose)
		if(returnierr == 0)
			cout << "Comm test successful." << endl;
		else
			cout << "Comm test failed." << endl;
  return(returnierr);
}

//======================================================================
template<typename PacketType, typename OrdinalType> 
int checkRankAndSize(Tpetra::Comm<PacketType, OrdinalType>& comm, int rank, int size, bool verbose, bool debug) {
  if(verbose) cout << "Checking rank and size... ";
  
  int ierr = 0;
  
  int comm_rank = comm.getMyImageID();
  if(comm_rank != rank)
    ierr++;
  int comm_size = comm.getNumImages();
  if(comm_size != size)
    ierr++;
  
  if(verbose)
    if(ierr == 0)
      cout << "Successful" << endl;
    else
      cout << "Failed" << endl;
  
  if(debug || (ierr != 0)) {
    cout << "rank = " << rank << ", myImageID = " << comm_rank << endl;
    cout << "size = " << size << ", numImages = " << comm_size << endl;
  }
  
  return(ierr);
}

//======================================================================
template<typename PacketType, typename OrdinalType>
int commMethods(Tpetra::Comm<PacketType, OrdinalType>& comm, bool verbose) {
	if(verbose) cout << "Starting Comm method testing..." << endl;
  
	OrdinalType count = Teuchos::OrdinalTraits<OrdinalType>::one();
	PacketType inVal;
	PacketType outVal;
  
  if(verbose) cout << "getMyImageID..." << endl;
  comm.getMyImageID(); // throw away output
  
  if(verbose) cout << "getNumImages..." << endl;
  comm.getNumImages(); // throw away output
  
  if(verbose) cout << "barrier..." << endl;
  comm.barrier();
  
	if(verbose) cout << "broadcast..." << endl;
	comm.broadcast(&inVal, count, 0);
  
	if(verbose) cout << "gatherAll..." << endl;
	comm.gatherAll(&inVal, &outVal, count);
  
	if(verbose) cout << "sumAll..." << endl;
	comm.sumAll(&inVal, &outVal, count);
  
	if(verbose) cout << "maxAll..." << endl;
	comm.maxAll(&inVal, &outVal, count);
  
	if(verbose) cout << "minAll..." << endl;
	comm.minAll(&inVal, &outVal, count);
  
	if(verbose) cout << "scanSum..." << endl;
	comm.scanSum(&inVal, &outVal, count);
  
  if(verbose) cout << "printInfo..." << endl;
  comm.barrier();
  if(verbose) comm.printInfo(cout); // only run if in verbose mode
  comm.barrier();
  
  if(verbose) cout << "Comm method testing finished." << endl;
	return(0);
}
