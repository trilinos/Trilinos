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

// Tpetra SparseMatrix tester

#include "Tpetra_ConfigDefs.hpp" // for <iostream> and <stdlib>
#include "Tpetra_SparseMatrix.hpp"
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

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

	// call test routine
	int ierr = 0;
	if(verbose) cout << "Starting SparseMatrixTest..." << endl;
	ierr += unitTests<int, float>(verbose, debug);
	ierr += unitTests<int, double>(verbose, debug);

	// finish up
	if(verbose)
		if(ierr == 0)
			cout << "SparseMatrix test successful." << endl;
		else
			cout << "SparseMatrix test failed." << endl;
	return(ierr);
}


//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
	int ierr = 0;
	int returnierr = 0;
	char const * OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
	char const * STName = Teuchos::ScalarTraits<ScalarType>::name();

	if(verbose) cout << "Starting unit tests for SparseMatrix<" << OTName << "," << STName << ">." << endl;

	// ======================================================================
	// code coverage section - just call functions, no testing
	// ======================================================================

	if(verbose) cout << "Constructors..." << endl;
	// default constructor
	Tpetra::SparseMatrix<OrdinalType, ScalarType> sm;
	// copy constructor
	Tpetra::SparseMatrix<OrdinalType, ScalarType> smClone(sm);

	if(verbose) cout << "Code coverage section finished." << endl;

	// ======================================================================
	// actual testing section - affects return code
	// ======================================================================

	if(verbose) cout << "Starting actual testing section... (none to do)" << endl;

	// finish up
	if(verbose)
		if(returnierr == 0)
			cout << "SparseMatrixTest <" << OTName << ", " << STName << "> passed." << endl;
		else
			cout << "SparseMatrixTest <" << OTName << ", " << STName << "> failed." << endl;
  
	return(returnierr);
}
