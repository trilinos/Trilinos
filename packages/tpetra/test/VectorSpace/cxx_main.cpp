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

// Tpetra VectorSpace tester
// Modified: 06-Feb-2003

#include <iostream>
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Version.hpp"

// function prototype
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug);

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

	if(verbose)
		cout << Tpetra::Tpetra_Version() << endl << endl;

	// call test routine
	int ierr = 0;
	if(verbose) cout << "Starting VectorSpaceTest..." << endl;
	ierr += unitTests<int, float>(verbose, debug);
	ierr += unitTests<int, double>(verbose, debug);

	// finish up
	if(verbose) 
		if(ierr == 0)
			cout << "VectorSpace test passed." << endl;
		else
			cout << "VectorSpace test failed." << endl;
	return(ierr);
}

template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
	Tpetra::SerialPlatform<OrdinalType, OrdinalType> platformE;
	Tpetra::SerialPlatform<OrdinalType, ScalarType> platformV;
	if(verbose) cout << "Starting unit tests for VectorSpace<" 
			 << Teuchos::OrdinalTraits<OrdinalType>::name() << "," 
			 << Teuchos::ScalarTraits<ScalarType>::name() << ">." << endl;

	//
	// code coverage section - just call functions, no testing
	//
	if(verbose) cout << "Starting code coverage section..." << endl;
	// constructors
	if(verbose) cout << "Constructors..." << endl;
	// taking an ElementSpace
	Tpetra::ElementSpace<OrdinalType> elementspace(10, 2, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> vectorspace(elementspace, platformV);
	// taking a BlockElementSpace
	const Tpetra::BlockElementSpace<OrdinalType> blockelementspace(elementspace, 3);
	Tpetra::VectorSpace<OrdinalType, ScalarType> blockvectorspace(blockelementspace, platformV);
	// cpy ctr
	Tpetra::VectorSpace<OrdinalType, ScalarType> v2(vectorspace);

	// print
	if(debug) {
		cout << "Overloaded << operator..." << endl;
		cout << vectorspace << endl;
	}

	// attribute access
	if(verbose) cout << "Attribute access methods..." << endl;
	OrdinalType temp = 0;
	temp = vectorspace.getNumGlobalEntries();
	temp = vectorspace.getNumMyEntries();
	temp = vectorspace.getIndexBase();
	temp = vectorspace.getMinLocalIndex();
	temp = vectorspace.getMaxLocalIndex();
	temp = vectorspace.getMinGlobalIndex();
	temp = vectorspace.getMaxGlobalIndex();
	temp = 0;
	temp = vectorspace.getGlobalIndex(temp);
	temp = vectorspace.getLocalIndex(temp);
  bool tempB = vectorspace.isMyLocalIndex(temp);
  tempB = vectorspace.isMyGlobalIndex(temp);

	// vector creation
	if(verbose) cout << "Vector creation..." << endl;
	Tpetra::Vector<OrdinalType, ScalarType>* vecptr = vectorspace.createVector();
	temp = vecptr->getNumMyEntries();
	delete vecptr;

	// accessors to other classes
	if(verbose) cout << "Class member accessors..." << endl;
	if(verbose) v2.platform().printInfo(cout);
	temp = v2.comm().getNumImages();

	if(verbose) cout << "Code coverage section finished." << endl;

	//
	// actual testing section - affects return code
	//

	if(verbose) cout << "Starting actual testing section..." << endl;

	// isCompatible
	if(verbose) cout << "Testing isCompatible... ";
	const Tpetra::ElementSpace<OrdinalType> compatibleES(10,0, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> compatibleVS(compatibleES, platformV);
	if(!vectorspace.isCompatible(compatibleVS))
		ierr++;
	const Tpetra::ElementSpace<OrdinalType> differentES(15, 2, platformE);
	Tpetra::VectorSpace<OrdinalType, ScalarType> differentVS(differentES, platformV);
	if(differentVS.isCompatible(vectorspace))
		ierr++;
	if(verbose)
		if(ierr == 0) 
			cout << "Passed" << endl;
		else
			cout << "Failed" << endl;
	returnierr += ierr;
	ierr = 0;

	// isSameAs
	if(verbose) cout << "Testing isSameAs... ";
	bool same = vectorspace.isSameAs(v2);
	if(!same)
		ierr++;
	same = vectorspace.isSameAs(differentVS);
	if(same)
		ierr++;
	if(verbose)
		if(ierr == 0) 
			cout << "Passed" << endl;
		else
			cout << "Failed" << endl;
	returnierr += ierr;
	ierr = 0;

  // assignment operator
  if(verbose) cout << "Testing asignment operator... ";
  same = vectorspace.isSameAs(differentVS);
  if(same)
    ierr++;
  differentVS = vectorspace;
  same = vectorspace.isSameAs(differentVS);
  if(!same)
    ierr++;
	if(verbose)
		if(ierr == 0) 
			cout << "Passed" << endl;
		else
			cout << "Failed" << endl;
	returnierr += ierr;
	ierr = 0;
	
	// finish up
	//char const * OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
	//char const * STName = Teuchos::ScalarTraits<ScalarType>::name();
	if(verbose)
		if(returnierr == 0)
			cout << "VectorSpaceTest <" 
			     << Teuchos::OrdinalTraits<OrdinalType>::name() << ", " 
			     << Teuchos::ScalarTraits<ScalarType>::name()  << "> passed." << endl;
		else
			cout << "VectorSpaceTest <" 
			     << Teuchos::OrdinalTraits<OrdinalType>::name() << ", " 
			     << Teuchos::ScalarTraits<ScalarType>::name() << ">failed." << endl;
	return(returnierr);
}
