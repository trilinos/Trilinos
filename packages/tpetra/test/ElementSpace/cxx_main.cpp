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
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_Version.hpp"

template<typename OrdinalType>
void esTester(bool verbose, bool debug, int size, int rank);
template<typename OrdinalType>
void isLgetG(OrdinalType low, OrdinalType high, Tpetra::ElementSpace<OrdinalType>& es);
template<typename OrdinalType>
void isGgetL(OrdinalType low, OrdinalType high, Tpetra::ElementSpace<OrdinalType>& es);

int main(int argc, char* argv[]) {
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
    << "Starting ElementSpaceTest..." << endl
    << Tpetra::Tpetra_Version() << endl
    << "****************************************\n";
  }
  
	esTester<int>(verbose, debug, size, rank);

#ifdef TPETRA_MPI
  MPI_Finalize();
#endif // TPETRA_MPI
  
	return(0);
}

template<typename OrdinalType>
void esTester(bool verbose, bool debug, int size, int rank) {
	OrdinalType nME = 5;
  OrdinalType nGE = nME * size;
  
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
#endif // TPETRA_MPI
	
	if(verbose) cout << "Creating es1(contiguous, tpetra-defined)...";
	Tpetra::ElementSpace<OrdinalType> es1(nGE, 2, platform);
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << es1 << endl;
	
	if(verbose) cout << "Creating es2(contiguous, user-defined)...";
	Tpetra::ElementSpace<OrdinalType> es2(-1, nME, 2, platform); // should be same as es1
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << es2 << endl;
	
	if(verbose) cout << "Testing isSameAs (contig)...";
	assert(es1.isSameAs(es2) == true);
	Tpetra::ElementSpace<OrdinalType> es2a((nGE * 2), (nME * 2), 3, platform); // should be different than es1 & es2
	assert(es1.isSameAs(es2a) == false);
	if(verbose) cout << "Successful." << endl;
	
	if(verbose) cout << "Creating es3(noncontiguous)...";
  // create list of GIDs; needs to be unique for each image
  // should be {0, 10, 20, ...} on node 0, {11, 21, 31} on node 1, etc.
  std::vector<OrdinalType> eList(nME); // allocate to size nME
  for(OrdinalType i = 0; i < nME; i++)
    eList[i] = 10 * (i + rank) + rank;
  
	Tpetra::ElementSpace<OrdinalType> es3(-1, nME, eList, 0, platform);
	if(verbose) cout << "Successful." << endl;
	if(debug) cout << es3 << endl;
	
	if(debug) {
		cout << "Testing isMyLID and getGID" << endl;
		isLgetG(0, 13, es3); 
		cout << "Testing isMyGID and getLID" << endl; 
		isGgetL(0, 60, es3);
	}
	
	if(verbose) cout << "Testing isSameAs (noncontig)...";
	Tpetra::ElementSpace<OrdinalType> es3a(nGE, nME, eList, es3.getIndexBase(), platform); // should be same as es3
	assert(es3.isSameAs(es3a) == true);
	eList[(nME / 2)] *= 10;
	Tpetra::ElementSpace<OrdinalType> es3b(nGE, nME, eList, es3.getIndexBase(), platform); // should be different than es3 & es3a
	assert(es3.isSameAs(es3b) == false);
	if(verbose) cout << "Successful." << endl;

	if(verbose) cout << "Testing copy constructor...";
	Tpetra::ElementSpace<OrdinalType> es4(es3);
	assert(es3.isSameAs(es4) == true);
	assert(es4.isSameAs(es3b) == false);
	if(verbose) cout << "Successful." << endl;

  if(verbose) cout << "Testing assignment operator...";
  assert(es3.isSameAs(es3b) == false);
  es3b = es3;
  assert(es3.isSameAs(es3b) == true);
  if(verbose) cout << "Successful." << endl;

#ifdef TPETRA_MPI
  if(verbose) cout << "Running in MPI mode, not testing getRemoteIDList." << endl;
#else
	if(verbose) cout << "Testing getRemoteIDList...";
	int const len = 4;
  std::vector<OrdinalType> gList(len);
  gList[0] = eList[1]; 
  gList[1] = eList[3]; 
  gList[2] = eList[5];
  gList[3] = eList[0];
  std::vector<OrdinalType> pList(len, 5);
  std::vector<OrdinalType> lList(len, 0);
	es3.getRemoteIDList(len, gList, pList, lList);
	if(debug) cout << "\nGID PID LID getLID" << endl;
	for(int i = 0; i < len; i++) {
		if(debug) cout << setw(3) << gList[i] << setw(4) << pList[i] << setw(4) << lList[i] << setw(4) << es3.getLID(gList[i]) << endl;
		assert(lList[i] == es3.getLID(gList[i]));
	}
	if(verbose) cout << "Successful." << endl;
#endif // TPETRA_MPI
	
	if(verbose) cout << "ElementSpace test successful." << endl;
}

template<typename OrdinalType>
void isLgetG(OrdinalType low, OrdinalType high, Tpetra::ElementSpace<OrdinalType>& es) {
	  for(OrdinalType i = low; i < high; i++) {
		  if(es.isMyLID(i))
			  cout << "LID" << setw(3) << i << " getGID? " << es.getGID(i) << endl;
	  }
}

template<typename OrdinalType>
void isGgetL(OrdinalType low, OrdinalType high, Tpetra::ElementSpace<OrdinalType>& es) {
	  for(OrdinalType i = low; i < high; i++) {
		  if(es.isMyGID(i))
			  cout << "GID" << setw(3) << i << " getLID? " << es.getLID(i) << endl;
	  }
}
