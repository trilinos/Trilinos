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

// Tpetra CisMatrix tester

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CombineMode.hpp"
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

// function prototype
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug);
template <typename OrdinalType, typename ScalarType>
int mapTester();

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

	// call test routine
	int ierr = 0;
	if(verbose) cout << "Starting CisMatrixTest..." << endl;
	//ierr += unitTests<int, float>(verbose, debug);
	ierr += unitTests<int, double>(verbose, debug);
  //int tmp = mapTester<int, double>();

	// finish up
	if(verbose)
		if(ierr == 0)
			cout << "CisMatrix test successful." << endl;
		else
			cout << "CisMatrix test failed." << endl;
	return(ierr);
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int mapTester() {
  typedef std::map<OrdinalType, ScalarType> OrdScalMap;
  typedef std::map<OrdinalType, OrdScalMap> MapOfMaps;

  OrdScalMap map1;
  map1[1] = 10.0;
  map1[2] = 20.0;
  map1[3] = 30.0;
  cout << "map1" << endl;
  for(typename OrdScalMap::iterator i = map1.begin(); i != map1.end(); i++)
    cout << "Key: " << (*i).first << " Data: " << (*i).second << endl;

  OrdScalMap map2;
  map2[1] = 11.0;
  map2[2] = 22.0;
  map2[3] = 33.0;
  cout << "map2" << endl;
  for(typename OrdScalMap::iterator i = map2.begin(); i != map2.end(); i++)
    cout << "Key: " << (*i).first << " Data: " << (*i).second << endl;

  MapOfMaps map3;
  map3[1] = map1;
  map3[2] = map2;
  cout << "map3" << endl;
  for(typename MapOfMaps::iterator i = map3.begin(); i != map3.end(); i++) {
    cout << "row/column: " << (*i).first << endl;
    OrdScalMap& innermap = (*i).second;
    for(typename OrdScalMap::iterator j = innermap.begin(); j != innermap.end(); j++)
      cout << "(" << (*j).first << "," << (*j).second << ")" << endl;
  }
}

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
	char const * OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
	char const * STName = Teuchos::ScalarTraits<ScalarType>::name();

	if(verbose) cout << "Starting unit tests for CisMatrix<" << OTName << "," << STName << ">." << endl;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	if(verbose) cout << "Constructors..." << endl;
	// have to create ElementSpace and VectorSpace first
  const Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;
	const Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
	Tpetra::ElementSpace<OrdinalType> elementspace(10, 0, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
  // constructor taking one VectorSpace
	Tpetra::CisMatrix<OrdinalType, ScalarType> sm(vectorspace);
  // constructor taking two VectorSpaces
  Tpetra::CisMatrix<OrdinalType, ScalarType> sm2(vectorspace, vectorspace);
	// copy constructor
	Tpetra::CisMatrix<OrdinalType, ScalarType> smClone(sm);

  // submissions
  if(verbose) cout << "Submitting entries..." << endl;
  sm.submitEntry(Tpetra::Replace, 0, 5, 1);
  sm.submitEntry(Tpetra::Replace, 0, 2, 2);
  sm.submitEntry(Tpetra::Replace, 1, 8, 0);
  sm.submitEntry(Tpetra::Replace, 1, 6, 3);
  sm.submitEntry(Tpetra::Replace, 2, 3, 2);
  sm.submitEntry(Tpetra::Replace, 3, 4, 0);
  sm.submitEntry(Tpetra::Replace, 3, 11, 1);
  sm.submitEntry(Tpetra::Replace, 3, 1, 2);
  sm.submitEntry(Tpetra::Add, 3, 1, 1);

  if(debug) cout << sm << endl;

  // scale
  if(verbose) cout << "scale..." << endl;
  sm.scale(2.0);
  if(debug) cout << sm << endl;

  // setAllToScalar
  if(verbose) cout << "setAllToScalar..." << endl;
  sm.setAllToScalar(6.0);
  if(debug) cout << sm << endl;

  // get attributes
  if(verbose) cout << "getNumGlobalNonzeros..." << endl;
  sm.getNumGlobalNonzeros(); // throw away output
  if(verbose) cout << "getNumMyNonzeros..." << endl;
  sm.getNumMyNonzeros(); // throw away output
  if(verbose) cout << "getNumGlobalDiagonals..." << endl;
  sm.getNumGlobalDiagonals(); // throw away output
  if(verbose) cout << "getNumMyDiagonals..." << endl;
  sm.getNumMyDiagonals(); // throw away output
  
  // retrieve VectorSpaces
  if(verbose) cout << "Retrieving VectorSpaces..." << endl;
  if(verbose) cout << "sm.getPrimaryDist()" << endl;
  sm.getPrimaryDist(); // throw away output
  //if(verbose) cout << "sm.getSecondaryDist()" << endl;
  //sm.getSecondaryDist(); // throw away output
  //if(verbose) cout << "sm.getDomainMap()" << endl;
  //sm.getDomainDist(); // throw away output
  //if(verbose) cout << "sm.getRangeMap()" << endl;
  //sm.getRangeDist(); // throw away output

  // fillComplete
  if(verbose) cout << "fillComplete..." << endl;
  sm.fillComplete();

  // print
  cout << sm << endl;
  

	if(verbose) cout << "Code coverage section finished." << endl;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	if(verbose) cout << "Starting actual testing section... (none to do)" << endl;

	// finish up
	if(verbose)
		if(returnierr == 0)
			cout << "CisMatrixTest <" << OTName << ", " << STName << "> passed." << endl;
		else
			cout << "CisMatrixTest <" << OTName << ", " << STName << "> failed." << endl;
  
	return(returnierr);
}
