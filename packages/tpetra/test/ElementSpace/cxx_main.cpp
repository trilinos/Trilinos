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
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif // TPETRA_MPI

template<typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

template<typename OrdinalType>
int esTest(bool verbose, bool debug, int myImageID, int numImages, OrdinalType const indexBase, bool const global, bool const contiguous, Tpetra::ElementSpace<OrdinalType> const& es);

template<typename OrdinalType>
int testLIDGID(Tpetra::ElementSpace<OrdinalType> const& es, OrdinalType nME, std::vector<OrdinalType> const& myGIDs);

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

  int myImageID = 0; // assume we are on serial
  int numImages = 1; // if MPI, will be reset later
  
  // initialize MPI if needed
#ifdef TPETRA_MPI
  numImages = -1;
  myImageID = -1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numImages);
  MPI_Comm_rank(MPI_COMM_WORLD, &myImageID);
  MPI_Barrier(MPI_COMM_WORLD);
#endif // TPETRA_MPI
  
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = (verbose && (myImageID == 0));
  
  // start the testing
	if(verbose) outputStartMessage("ElementSpace");
  int ierr = 0;
  
  //mpiBreakpoint(myImageID);

  // call the actual test routines
	ierr += unitTests<int>(verbose, debug, myImageID, numImages);
  
	// finish up
#ifdef TPETRA_MPI
  MPI_Finalize();
#endif
	if(verbose) outputEndMessage("ElementSpace", (ierr == 0));
	return(ierr);
}

//======================================================================
template<typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  std::string className = "ElementSpace<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  int ierr = 0;
  int returnierr = 0;

  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();
  OrdinalType const negOne = zero - one;
  
  // ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
  Tpetra::MpiComm<OrdinalType, OrdinalType> comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
  Tpetra::SerialComm<OrdinalType, OrdinalType> comm;
#endif // TPETRA_MPI

  // constructor #1 - nGE, iB, Platform
  if(verbose) cout << "Calling ElementSpace constructor #1 (contiguous, tpetra-defined)..." << endl;
  comm.barrier();
  Tpetra::ElementSpace<OrdinalType> es1(intToOrdinal<OrdinalType>(5*numImages), zero, platform);
  
  // constructor #2 - nGE, nME, iB, Platform
  if(verbose) cout << "Calling ElementSpace constructor #2 (contiguous, user-defined)..." << endl;
  comm.barrier();
  Tpetra::ElementSpace<OrdinalType> es2(negOne, intToOrdinal<OrdinalType>(5), zero, platform);
  
	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

  // fixtures
  OrdinalType const iB = zero; // magic number for indexBase, but normal value
	OrdinalType const nME = intToOrdinal<OrdinalType>(5); // magic number for numMyElements
  OrdinalType const nGE = nME * intToOrdinal<OrdinalType>(numImages);
  std::vector<OrdinalType> myGIDs;
  std::vector<OrdinalType> tempVector;
  std::vector<OrdinalType> expected;

  generateColumn(myGIDs, myImageID, nME); // give each image a column from the generator (length nME)
  
  // ========================================
  // test ctr #3 - non-contig
  // ========================================
  if(verbose) cout << "Testing ElementSpace constructor #3 (noncontiguous)...";
  comm.barrier();
  Tpetra::ElementSpace<OrdinalType> elementspace(nGE, nME, myGIDs, iB, platform);
  ierr = ctrTest(elementspace, myGIDs, nGE, generateAreaMin<OrdinalType>(0,0,numImages-1,nME-1),
                 generateAreaMax<OrdinalType>(0,0,numImages-1,nME-1), iB, (numImages > 1), false); // should be global if we're running on more than one image
  if(ierr != 0) {
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  // call isMyLID, isMyGID, getGID, and getLID for every element we own
  if(verbose) cout << "Testing isMyLID, isMyGID, getLID, and getGID... ";
  ierr = testLIDGID(elementspace, nME, myGIDs);
  if(ierr != 0) {
    if(verbose) cout << "failed" << endl;
    ierr = 1;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  /*
  // test getRemoteIDList
  if(verbose) cout << "Testing getRemoteIDList... ";
  std::vector<OrdinalType> allGIDs(nGE);
  generateMultipleColumns(allGIDs, 0, numImages-1, nME);
  std::vector<OrdinalType> imageIDList(nGE, intToOrdinal<OrdinalType>(-99));
  std::vector<OrdinalType> LIDList(nGE, intToOrdinal<OrdinalType>(-99));
	elementspace.getRemoteIDList(allGIDs, imageIDList, LIDList);

  generateXCoords(allGIDs, expected);
  generateYCoords(allGIDs, tempVector);
  if(debug) {
    if(verbose) cout << endl;
    outputData(myImageID, numImages, "allGIDs: " + Tpetra::toString(allGIDs));
    outputData(myImageID, numImages, "imageIDList: " + Tpetra::toString(imageIDList));
    outputData(myImageID, numImages, "expected:    " + Tpetra::toString(expected));
    outputData(myImageID, numImages, "LIDList:  " + Tpetra::toString(LIDList));
    outputData(myImageID, numImages, "expected: " + Tpetra::toString(tempVector));
    //cout << elementspace << endl;
    if(verbose) cout << "getRemoteIDList test ";
  }
  if(imageIDList != expected || LIDList != tempVector) {
    ierr++;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;
  */
  
  // test isSameAs (contig)
  comm.barrier();
	if(verbose) cout << "Testing isSameAs (contig)... ";
	Tpetra::ElementSpace<OrdinalType> es2a(nGE, nME, one, platform); // should be different than es1 & es2
  /*if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << es1 << endl;
    comm.barrier();
    cout << es2 << endl;
    comm.barrier();
    cout << es2a << endl;
    comm.barrier();
    if(verbose) cout << "isSameAs (contig) test ";
    }*/
  if(!es1.isSameAs(es2) || es2.isSameAs(es2a)) { // es1 and es2 should be the same, es2 and es2a should be different
    ierr++;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;


  // test isSameAs (noncontig)
	if(verbose) cout << "Testing isSameAs (noncontig)... ";
	Tpetra::ElementSpace<OrdinalType> es3a(nGE, nME, myGIDs, elementspace.getIndexBase(), platform); // should be same as es3
  tempVector = myGIDs;
	tempVector[(nME / intToOrdinal<OrdinalType>(2))] *= intToOrdinal<OrdinalType>(10); // change one of the GID values given as myElementList
	Tpetra::ElementSpace<OrdinalType> es3b(nGE, nME, tempVector, elementspace.getIndexBase(), platform); // should be different than es3 & es3a
  /*if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << es3a << endl;
    comm.barrier();
    cout << es3b << endl;
    comm.barrier();
    if(verbose) cout << "isSameAs (noncontig) test ";
    }*/
  if(!elementspace.isSameAs(es3a) || elementspace.isSameAs(es3b)) { // es3 and es3a should be the same, es3 and es3b should be different
    ierr++;
    if(verbose) cout << "failed" << endl;
  }
  else
    if(verbose) cout << "passed" << endl;
  returnierr += ierr;
  ierr = 0;

  
  // test cpy ctr
	if(verbose) cout << "Testing copy constructor... ";
	Tpetra::ElementSpace<OrdinalType> esClone(elementspace);
  /*if(debug) {
    if(verbose) cout << endl;
    comm.barrier(); 
    cout << elementspace << endl;
    comm.barrier();
    cout << esClone << endl; 
    comm.barrier();
    if(verbose) cout << "Copy Constructor test ";
  }*/
	if(!elementspace.isSameAs(esClone) || esClone.isSameAs(es3b)) {
    if(verbose) cout << "failed" << endl;
    ierr++;
  }
  else
    if(verbose) cout << "passed" << endl;
	returnierr += ierr;
  ierr = 0;


  // test operator=
  if(verbose) cout << "Testing assignment operator... ";
  assert(elementspace.isSameAs(es3b) == false); // if this fails, the problem is in this tester, not in ElementSpace.
  es3b = elementspace;
  /*if(debug) {
    if(verbose) cout << endl;
    comm.barrier();
    cout << elementspace << endl;
    comm.barrier();
    cout << es3b << endl;
    comm.barrier();
    if(verbose) cout << "Operator = test ";
  }*/
  if(!elementspace.isSameAs(es3b)) {
    ierr++;
    if(verbose) cout << "failed" << endl;
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


// test properties of ElementSpace. This is easy to do, but tedious
// so it's done here
//======================================================================
template<typename OrdinalType>
int ctrTest(Tpetra::ElementSpace<OrdinalType> const& es, std::vector<OrdinalType> const& myGIDs, 
            OrdinalType const numGlobalElements, OrdinalType const minAllGID, OrdinalType const maxAllGID, 
            OrdinalType const indexBase, bool const global, bool const contiguous)
{
  int ierr = 0;

  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const one = Teuchos::OrdinalTraits<OrdinalType>::one();

  // compute nME, maxMyGID, maxAllGID, minMyGID, minAllGID
  OrdinalType numMyElements = myGIDs.size();
  OrdinalType minMyGID = *min_element(myGIDs.begin(), myGIDs.end());
  OrdinalType maxMyGID = *max_element(myGIDs.begin(), myGIDs.end());

  OrdinalType esOrd; // temp variable for functions returning an OrdinalType
  bool esBool; // temp variable for functions returning a boolean
  
  esOrd = es.getNumGlobalElements();
  if(esOrd != numGlobalElements) 
    ierr++;

  esOrd = es.getNumMyElements();
  if(esOrd != numMyElements)
    ierr++;

  esOrd = es.getIndexBase();
  if(esOrd != indexBase)
    ierr++;

  esOrd = es.getMinLID();
  if(esOrd != zero)
    ierr++;

  esOrd = es.getMaxLID();
  if(esOrd != (numMyElements - one))
    ierr++;

  esOrd = es.getMinMyGID();
  if(esOrd != minMyGID)
    ierr++;

  esOrd = es.getMaxMyGID();
  if(esOrd != maxMyGID)
    ierr++;

  esOrd = es.getMinAllGID();
  if(esOrd != minAllGID)
    ierr++;

  esOrd = es.getMaxAllGID();
  if(esOrd != maxAllGID)
    ierr++;

  esBool = es.isGlobal();
  if(esBool != global)
    ierr++;

  esBool = es.isContiguous();
  if(esBool != contiguous)
    ierr++;

  std::vector<OrdinalType> const& myGlobalElements = es.getMyGlobalElements();
  if(myGlobalElements != myGIDs)
    ierr++;
  
  return(ierr);
}

//======================================================================
template<typename OrdinalType>
int testLIDGID(Tpetra::ElementSpace<OrdinalType> const& es, OrdinalType nME, std::vector<OrdinalType> const& myGIDs) {
  OrdinalType const zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
  OrdinalType const negOne = zero - Teuchos::OrdinalTraits<OrdinalType>::one();
  int ierr = 0;

  for(OrdinalType i = zero; i < nME; i++) {
    if(!es.isMyLID(i))
      ierr++;
    if(es.getGID(i) != myGIDs[i])
      ierr++;
    if(!es.isMyGID(myGIDs[i]))
      ierr++;
    if(es.getLID(myGIDs[i]) != i)
      ierr++;
  }
  // call the border cases for isMyLID
  if(es.isMyLID(negOne))
    ierr++;
  if(es.isMyLID(nME))
    ierr++;

  return(ierr);
}
