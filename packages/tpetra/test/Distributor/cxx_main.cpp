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

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_Version.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiDistributor.hpp"
#else
#include "Tpetra_SerialDistributor.hpp"
#endif // TPETRA_MPI

// function prototype
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
    cout << "\n****************************************\n" 
    << "Starting DistributorTest..." << endl
    << Tpetra::Tpetra_Version() << endl
    << "****************************************\n";
  }
  int ierr = 0;
  
	ierr += unitTests<int, float>(verbose, debug, rank, size);
	ierr += unitTests<int, double>(verbose, debug, rank, size);

	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
  
	if(verbose) 
		if(ierr == 0)
			cout << "Distributor test successful." << endl;
		else
			cout << "Distributor test failed." << endl;
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int rank, int size) {
	if(verbose) cout << "Starting unit tests for Distributor<" 
									 << Teuchos::OrdinalTraits<OrdinalType>::name() << "," 
									 << Teuchos::ScalarTraits<ScalarType>::name() << ">." << endl;
	int ierr = 0;
	int returnierr = 0;
	
	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
	if(verbose) cout << "Starting code coverage section..." << endl;
	
	// constructors
	if(debug) cout << "Constructors..." << endl;
#ifdef TPETRA_MPI
  Tpetra::MpiDistributor<OrdinalType, OrdinalType> distrib(MPI_COMM_WORLD);
#else
  Tpetra::SerialDistributor<OrdinalType, OrdinalType> distrib;
#endif // TPETRA_MPI
  // cpy ctr
	
	// print
	if(debug) {
		cout << "Overloaded << operator..." << endl;
		cout << distrib << endl;
	}
	
  // assignment operator
  if(debug) cout << "assignment operator..." << endl;
  //v2 = vector;
	
	if(verbose) cout << "Code coverage section finished." << endl;
	
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================
	
	if(verbose) cout << "Starting actual testing section...(none to do)" << endl;
	
	// finish up
	if(verbose) {
    cout << "DistributorTest <" 
         << Teuchos::OrdinalTraits<OrdinalType>::name() << ", " 
         << Teuchos::ScalarTraits<ScalarType>::name() << "> ";
		if(returnierr == 0)
			cout << "passed." << endl;
		else
			cout << "failed." << endl;
  }
	return(returnierr);
}
