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
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_Import.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size);
template <typename OrdinalType>
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
	if(verbose) outputStartMessage("ImportExport");
  int ierr = 0;
  
  // call the actual test routines
	ierr += unitTests<int>(verbose, debug, rank, size);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("ImportExport", (ierr == 0));
	return(ierr);
}

//======================================================================
template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size) {
  std::string className = "Import<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  codeCoverage<OrdinalType>((verbose && debug), rank, size);
	
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  if(verbose && debug) outputSubHeading("Starting actual testing section... ");

#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
  Tpetra::MpiComm<OrdinalType, OrdinalType> comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
  Tpetra::SerialComm<OrdinalType, OrdinalType> comm;
#endif

  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const negOne = zero - Teuchos::OrdinalTraits<OrdinalType>::one();
  Tpetra::ElementSpace<OrdinalType> source(negOne, intToOrdinal<OrdinalType>(10), zero, platform);
  Tpetra::ElementSpace<OrdinalType> target(negOne, intToOrdinal<OrdinalType>(5), zero, platform);
  Tpetra::Import<OrdinalType> importer(source, target);

  // test that numSame + numPermute + numRemote = target.getNumMyElements()
  if(verbose) cout << "Testing same/permute/remote sum... ";

  OrdinalType same = importer.getNumSameIDs();
  OrdinalType permute = importer.getNumPermuteIDs();
  OrdinalType remote = importer.getNumRemoteIDs();
  OrdinalType sum = same + permute + remote;
  OrdinalType expectedSum = target.getNumMyElements();
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << "[Image " << rank << "] NumSameIDs:    " << same << endl;
    cout << "[Image " << rank << "] NumPermuteIDs: " << permute << endl;
    cout << "[Image " << rank << "] NumRemoteIDs:  " << remote << endl;
    cout << "[Image " << rank << "] Expected Sum:  " << expectedSum << endl;
    comm.barrier();
    if(verbose) cout << "same/permute/remote sum test ";
  }
  if(sum != expectedSum) {
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
template <typename OrdinalType>
void codeCoverage(bool verbose, int rank, int size) { 
  if(verbose) outputSubHeading("Starting code coverage section...");

  if(verbose) cout << "Creating Platform and Comm..." << endl;
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
  Tpetra::MpiComm<OrdinalType, OrdinalType> comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
  Tpetra::SerialComm<OrdinalType, OrdinalType> comm;
#endif

  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const negOne = zero - Teuchos::OrdinalTraits<OrdinalType>::one();
  OrdinalType const five = intToOrdinal<OrdinalType>(5);
  OrdinalType const ten = intToOrdinal<OrdinalType>(10);

  if(verbose) cout << "Creating ElementSpaces..." << endl;
  Tpetra::ElementSpace<OrdinalType> source(negOne, ten, zero, platform);
  Tpetra::ElementSpace<OrdinalType> target(negOne, five, zero, platform);

  if(verbose) cout << "Import constructor..." << endl;
  Tpetra::Import<OrdinalType> importer(source, target);

  if(verbose) cout << "Import copy constructor..." << endl;
  Tpetra::Import<OrdinalType> importer2(importer);

  if(verbose) cout << "getNumSameIDs..." << endl;
  importer.getNumSameIDs();

  if(verbose) cout << "getNumPermuteIDs..." << endl;
  importer.getNumPermuteIDs();

  if(verbose) cout << "getPermuteFromLIDs..." << endl;
  importer.getPermuteFromLIDs();

  if(verbose) cout << "getPermuteToLIDs..." << endl;
  importer.getPermuteToLIDs();

  if(verbose) cout << "getNumRemoteIDs..." << endl;
  importer.getNumRemoteIDs();

  if(verbose) cout << "getRemoteIDs..." << endl;
  importer.getRemoteLIDs();

  if(verbose) cout << "getNumExportIDs..." << endl;
  importer.getNumExportIDs();

  if(verbose) cout << "getExportLIDs..." << endl;
  importer.getExportLIDs();

  if(verbose) cout << "getExportImageIDs..." << endl;
  importer.getExportImageIDs();

  if(verbose) cout << "getSourceSpace..." << endl;
  importer.getSourceSpace();

  if(verbose) cout << "getTargetSpace..." << endl;
  importer.getTargetSpace();

  if(verbose) cout << "assignment operator..." << endl;
  importer2 = importer;
}

