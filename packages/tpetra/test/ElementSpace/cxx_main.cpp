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
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif // TPETRA_MPI

template<typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size);

template<typename OrdinalType>
void isLgetG(OrdinalType low, OrdinalType high, Tpetra::ElementSpace<OrdinalType>& es);
template<typename OrdinalType>
void isGgetL(OrdinalType low, OrdinalType high, Tpetra::ElementSpace<OrdinalType>& es);

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
  MPI_Barrier(MPI_COMM_WORLD);
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = (verbose && (rank == 0));
  
  // start the testing
	if(verbose) outputStartMessage("ElementSpace");
  int ierr = 0;
  
  // call the actual test routines
	ierr += unitTests<int>(verbose, debug, rank, size);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("ElementSpace", (ierr == 0));
	return(ierr);
}

//======================================================================
template<typename OrdinalType>
int unitTests(bool verbose, bool debug, int rank, int size) {
  std::string className = "ElementSpace<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;

	OrdinalType nME = intToOrdinal<OrdinalType>(5); // magic number for numMyElements
  OrdinalType nGE = intToOrdinal<OrdinalType>(5 * size); // magic number for numGlobalElements
  OrdinalType const zero = intToOrdinal<OrdinalType>(0); // magic number for indexBase
  OrdinalType const negOne = intToOrdinal<OrdinalType>(-1);
  
  // ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
  if(verbose) outputSubHeading("Starting code coverage section...");
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
  Tpetra::MpiComm<OrdinalType, OrdinalType> comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
  Tpetra::SerialComm<OrdinalType, OrdinalType> comm;
#endif // TPETRA_MPI
	// constructor #1 - nGE, iB, Platform
	if(verbose) cout << "Creating es1(contiguous, tpetra-defined)..." << endl;
	Tpetra::ElementSpace<OrdinalType> es1(nGE, zero, platform);
	if(debug) {comm.barrier(); cout << es1 << endl; comm.barrier();}
	// constructor #2 - nGE, nME, iB, Platform
	if(verbose) cout << "Creating es2(contiguous, user-defined)..." << endl;
	Tpetra::ElementSpace<OrdinalType> es2(negOne, nME, zero, platform); // should be same as es1
	if(debug) {comm.barrier(); cout << es2 << endl; comm.barrier();}
	// constructor #3 - nGE, nME, eList, iB, Platform
	if(verbose) cout << "Creating es3(noncontiguous)..." << endl;
  std::vector<OrdinalType> eList;   // give each image nME from generator
  generateColumn(eList, rank, nME);
	Tpetra::ElementSpace<OrdinalType> es3(negOne, nME, eList, zero, platform);
	if(debug) {comm.barrier(); cout << es3 << endl; comm.barrier();}

  // test isMyLID/isMyGID/getGID/getLID
	/*if(debug) {
		if(verbose) cout << "Testing isMyLID and getGID" << endl;
		isLgetG(0, 13, es3); 
		if(verbose) cout << "Testing isMyGID and getLID" << endl; 
		isGgetL(0, 60, es3);
    }*/

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  if(verbose) outputSubHeading("Starting actual testing section...");

  // test isSameAs (contig)
  OrdinalType const one = intToOrdinal<OrdinalType>(1); // magic number for new indexBase
	if(verbose) cout << "Testing isSameAs (contig)... ";
	Tpetra::ElementSpace<OrdinalType> es2a(nGE, nME, one, platform); // should be different than es1 & es2
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << es2a << endl;
    comm.barrier();
    if(verbose) cout << "isSameAs (contig) test ";
  }
  if(!es1.isSameAs(es2) || es2.isSameAs(es2a)) {
    ierr++;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // test isSameAs (noncontig)
	if(verbose) cout << "Testing isSameAs (noncontig)... ";
	Tpetra::ElementSpace<OrdinalType> es3a(nGE, nME, eList, es3.getIndexBase(), platform); // should be same as es3
	eList[(nME / intToOrdinal<OrdinalType>(2))] *= intToOrdinal<OrdinalType>(10); // change one of the GID values given as myElementList
	Tpetra::ElementSpace<OrdinalType> es3b(nGE, nME, eList, es3.getIndexBase(), platform); // should be different than es3 & es3a
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << es3a << endl;
    comm.barrier();
    cout << es3b << endl;
    comm.barrier();
    if(verbose) cout << "isSameAs (noncontig) test ";
  }
  if(!es3.isSameAs(es3a) || es3.isSameAs(es3b)) {
    ierr++;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // test cpy ctr
	if(verbose) cout << "Testing copy constructor... ";
	Tpetra::ElementSpace<OrdinalType> es4(es3);
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier(); 
    cout << es4 << endl; 
    comm.barrier();
    if(verbose) cout << "Copy Constructor test ";
  }
	if(!es3.isSameAs(es4) || es4.isSameAs(es3b)) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
	returnierr += ierr;
  ierr = 0;

  // test operator=
  if(verbose) cout << "Testing assignment operator... ";
  assert(es3.isSameAs(es3b) == false); // if this fails, the problem is in this tester, not in ElementSpace.
  es3b = es3;
  if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << es3b << endl;
    comm.barrier();
    if(verbose) cout << "Operator = test ";
  }
  if(!es3.isSameAs(es3b)) {
    ierr++;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // test getRemoteIDList
#ifdef TPETRA_MPI
  if(verbose) cout << "Running in MPI mode, not testing getRemoteIDList." << endl;
#else
	if(verbose) cout << "Testing getRemoteIDList... ";
  OrdinalType const len = intToOrdinal<OrdinalType>(4); // magic number for length of these arrays
  std::vector<OrdinalType> gList(len); 
  gList[0] = eList[1]; 
  gList[1] = eList[3]; 
  gList[2] = eList[5];
  gList[3] = eList[0];
  std::vector<OrdinalType> pList(len, intToOrdinal<OrdinalType>(5));
  std::vector<OrdinalType> lList(len, intToOrdinal<OrdinalType>(0));
	es3.getRemoteIDList(len, gList, pList, lList);
	if(debug) cout << "\nGID PID LID getLID" << endl;
	for(OrdinalType i = zero; i < len; i++) {
		if(debug) cout << setw(3) << gList[i] << setw(4) << pList[i] << setw(4) << lList[i] << setw(4) << es3.getLID(gList[i]) << endl;
		if(lList[i] != es3.getLID(gList[i]))
      ierr++;
	}
  if(debug) cout << "getRemoteIDList test ";
  if(ierr != 0) {
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
#endif // TPETRA_MPI
	
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
// functions for testing LIDs & GIDs
//======================================================================
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
