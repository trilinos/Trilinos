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
#include "Tpetra_Import.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_Version.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#include "Tpetra_MpiPlatform.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#endif

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

	if (verbose)
		cout << Tpetra::Tpetra_Version() << endl << endl;
	// call test routine
	int ierr = 0;
	if(verbose) cout << "Starting ImportExportTest..." << endl;
	ierr += unitTests<int, double>(verbose, debug);

	// finish up
	if(verbose) {
		if(ierr == 0)
			cout << "ImportExport test successful." << endl;
		else
			cout << "ImportExport test failed." << endl;
  }

	return(ierr);
};

//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug) {
  std::string OTName = Teuchos::OrdinalTraits<OrdinalType>::name();
  std::string STName = Teuchos::ScalarTraits<ScalarType>::name();

	if(verbose) cout << "Starting code coverage for Import<" << OTName << "," << STName << ">." << endl;

	// have to create ElementSpaces first
#ifdef TPETRA_MPI
  Tpetra::MpiPlatform<OrdinalType, OrdinalType> platform(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType, OrdinalType> platform;
#endif

	Tpetra::ElementSpace<OrdinalType> source(10, 0, platform);
  Tpetra::ElementSpace<OrdinalType> target(10, 0, platform);

  Tpetra::Import<OrdinalType> importer(source, target);
  if(debug) cout << importer << endl;

	if(verbose) cout << "Code coverage <" << OTName << ", " << STName << "> section finished." << endl;

	return(0);
}
