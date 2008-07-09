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
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_ConfigDefs.hpp" // for map, vector, string, and iostream 
// #include "Tpetra_Vector.hpp"
#ifdef HAVE_MPI
#include "Tpetra_MpiPlatform.hpp"
#else
# include "Tpetra_SerialPlatform.hpp"
#endif // HAVE_MPI
#include "Tpetra_Map.hpp"

using namespace Teuchos;

template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages);

#define TEST_C1_WITH_INVALID_ARGS(nG,iB,str)          \
{                                                     \
    bool caught_expected_error = false;               \
    try {                                             \
      Tpetra::Map<OrdinalType> bad_map(nG,iB,platform);      \
    }                                                 \
    catch (std::invalid_argument &ia) {               \
      caught_expected_error = true;                   \
      if(verbose){                                    \
        cout << str << " threw expected error:" << endl; \
      }                                               \
      if(debug){                                      \
        cout << ia.what() << endl;                    \
      }                                               \
    }                                                 \
    if (!caught_expected_error){                      \
      ierr++;                                         \
      if(verbose){                                    \
         cout << str << " DID NOT catch expected error." << endl; \
      }                                               \
    }                                                 \
}

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
  if(verbose) outputStartMessage("Map");
  int ierr = 0;

  // call the actual test routines
  ierr += unitTests<int, int>(verbose, debug, myImageID, numImages);
  ierr += unitTests<int, float>(verbose, debug, myImageID, numImages);
  ierr += unitTests<int, double>(verbose, debug, myImageID, numImages);

  // finish up
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  if(verbose) outputEndMessage("Map", (ierr == 0));
  return(ierr);
}


//======================================================================
template <typename OrdinalType, typename ScalarType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  std::string className = "Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ", " + Teuchos::ScalarTraits<ScalarType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  (void)debug;
  (void)myImageID;
  (void)numImages;

  int returnierr = 0;

  // ======================================================================
  // code coverage section - just call functions, no testing
  // ======================================================================
  
  // create Platform and Comm
#ifdef HAVE_MPI
  Tpetra::MpiPlatform<OrdinalType> platform(MPI_COMM_WORLD);
#else
  Tpetra::SerialPlatform<OrdinalType> platform;
#endif
  RCP< Comm<OrdinalType> > comm = platform.createComm();

  int ierr = 0;

  // test error catching
  {
    // Constructor 1
    // numGlobalEntries must be >= 0 and consistent across all procs
    // indexBase must be consistent across all procs
    if(verbose) cout << "Invalid arguments to constructor (numGlobal,indexBase,platform)..." << endl;
    TEST_C1_WITH_INVALID_ARGS(-1,0,"Negative numGlobalEntries...");
#ifdef HAVE_MPI
    TEST_C1_WITH_INVALID_ARGS((myImageID == 0 ? 0 : 1),0,"Inconsistent numGlobalEntries...");
    TEST_C1_WITH_INVALID_ARGS(0,(myImageID == 0 ? 0 : 1),"Inconsistent numGlobalEntries...");
#endif
  }

  // ======================================================================
  // actual testing section - affects return code
  // ======================================================================
  {

  }

  if(verbose) {
    if(returnierr == 0) {
      outputHeading("Unit tests for " + className + " passed.");
    }
    else {
      outputHeading("Unit tests for " + className + " failed.");
    }
  }
  return(returnierr);
}
