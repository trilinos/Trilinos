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
#include "Tpetra_Version.hpp"
#include "Tpetra_SerialPlatform.hpp"

#ifdef TPETRA_MPI
#include "Tpetra_MpiPlatform.hpp"
template <typename OrdinalType, typename ScalarType>
int mpiTests(bool verbose, bool debug);
#endif // TPETRA_MPI

template <typename OrdinalType, typename ScalarType>
int serialTests(bool verbose, bool debug);

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
    cout << "\n****************************************\n" 
    << "Starting PlatformTest..." << endl
    << Tpetra::Tpetra_Version() << endl
    << "****************************************\n";
  }
  int ierr = 0;
  
  // test SerialPlatform
	ierr += serialTests<int, double>(verbose, debug);
  
  // test MpiPlatform (if enabled)
#ifdef TPETRA_MPI
  if(verbose) cout << "Testing MPI functionality..." << endl;
	ierr += mpiTests<int, double>(verbose, debug);
  MPI_Finalize();
#endif
  
	// finish up
	if(verbose) 
		if(ierr == 0)
			cout << "Platform test passed." << endl;
		else
			cout << "Platform test failed." << endl;
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int serialTests(bool verbose, bool debug) {
	if(verbose) cout << "Creating SerialPlatform object...";
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platform;
	if(debug) cout << platform << endl;
	if(verbose) cout << "Successful." << endl;
	
	if(verbose) cout << "Creating SerialComm objects...";
	Tpetra::Comm<ScalarType, OrdinalType>* comm1 = platform.createScalarComm();
	Tpetra::Comm<OrdinalType, OrdinalType>* comm2 = platform.createOrdinalComm();
	delete comm1;
	delete comm2;
	if(verbose) cout << "Successful." << endl;
	
	if(verbose) cout << "Creating SerialDistributor objects...";
	Tpetra::Distributor<ScalarType, OrdinalType>* distributor1 = platform.createScalarDistributor();
	Tpetra::Distributor<OrdinalType, OrdinalType>* distributor2 = platform.createOrdinalDistributor();
	delete distributor1;
	delete distributor2;
	if(verbose) cout << "Successful." << endl;
  
  return(0);
}

//======================================================================
#ifdef TPETRA_MPI
template <typename OrdinalType, typename ScalarType>
int mpiTests(bool verbose, bool debug) {
	if(verbose) cout << "Creating MpiPlatform object...";
	Tpetra::MpiPlatform<OrdinalType, ScalarType> platform2(MPI_COMM_WORLD);
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Creating MpiComm objects...";
	Tpetra::Comm<ScalarType, OrdinalType>* comm3 = platform2.createScalarComm(); 
	Tpetra::Comm<OrdinalType, OrdinalType>* comm4 = platform2.createOrdinalComm();
  delete comm3;
  delete comm4;
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Creating MpiDistributor objects...";
	Tpetra::Distributor<ScalarType, OrdinalType>* distributor3 = platform2.createScalarDistributor();
	Tpetra::Distributor<OrdinalType, OrdinalType>* distributor4 = platform2.createOrdinalDistributor();
	delete distributor3;
	delete distributor4;
	if(verbose) cout << "Successful." << endl;

	return(0);
}
#endif // TPETRA_MPI
