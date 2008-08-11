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
#include <Teuchos_RCP.hpp>
#ifdef HAVE_MPI
# include "Tpetra_MpiPlatform.hpp"
#else
# include "Tpetra_SerialPlatform.hpp"
#endif // HAVE_MPI
#include <Teuchos_CommandLineProcessor.hpp>

using namespace Teuchos;

template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

int main(int argc, char* argv[]) {
  int myImageID = 0; // assume we are on serial
  int numImages = 1; // if MPI, will be reset later

  // initialize MPI if needed
#ifdef HAVE_MPI
  numImages = -1;
  myImageID = -1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numImages);
  MPI_Comm_rank(MPI_COMM_WORLD, &myImageID);
#endif // HAVE_MPI

  bool verbose = false; 
  bool debug = false;
  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet"  ,&verbose,"Print messages and results.");
  cmdp.setOption("debug"  ,"nodebug",&debug  ,"Print debugging info.");
  cmdp.parse(argc,argv);
  // change verbose to only be true on Image 0
  // if debug is enabled, it will still output on all nodes
  verbose = ((verbose || debug) && (myImageID == 0));

  // start the testing
  if(verbose) outputStartMessage("Platform");
  int ierr = 0;

  // call the actual test routines
  ierr += unitTests<int>(verbose, debug, myImageID, numImages);

  // finish up
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  if(verbose) outputEndMessage("Platform", (ierr == 0));
  return(ierr);
}

//======================================================================
template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  (void)debug; 
  (void)myImageID;
  (void)numImages;
  std::string className = "Platform<" + OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  //int ierr = 0; // unused now that there's no actual tests
  int returnierr = 0;

  // ======================================================================
  // code coverage section - just call functions, no testing
  // ======================================================================

#ifdef HAVE_MPI
  // default constructor
  if(verbose) cout << "MpiPlatform default constructor..." << endl;
  Tpetra::MpiPlatform<OrdinalType> platform( rcp(new OpaqueWrapper<MPI_Comm>(MPI_COMM_WORLD)) );
  // copy constructor
  if(verbose) cout << "MpiPlatform copy constructor..." << endl;
  Tpetra::MpiPlatform<OrdinalType> platform2(platform);
#else
  // default constructor
  if(verbose) cout << "SerialPlatform default constructor..." << endl;
  Tpetra::SerialPlatform<OrdinalType> platform;
  // copy constructor
  if(verbose) cout << "SerialPlatform copy constructor..." << endl;
  Tpetra::SerialPlatform<OrdinalType> platform2(platform);
#endif

  // clone
  if(verbose) cout << "clone..." << endl;
  RCP< Tpetra::Platform<OrdinalType> > platform3 = platform.clone();

  // createComm
  if(verbose) cout << "createComm..." << endl;
  RCP< Comm<OrdinalType> > comm1 = platform.createComm();

  // ======================================================================
  // actual testing section - affects return code
  // ======================================================================

  // ... none to do ...

  // ======================================================================
  // finish up
  // ======================================================================

  comm1->barrier();
  if(verbose) {
    if(returnierr == 0)
      outputHeading("Unit tests for " + className + " passed.");
    else
      outputHeading("Unit tests for " + className + " failed.");
  }
  return(returnierr);
}
