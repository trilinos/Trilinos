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
#include <Teuchos_FancyOStream.hpp>

using namespace Teuchos;

template <typename OrdinalType>
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
        *fos << ia.what() << endl;                    \
      }                                               \
    }                                                 \
    if (!caught_expected_error){                      \
      ierr++;                                         \
      if(verbose){                                    \
         cout << str << " DID NOT catch expected error." << endl; \
      }                                               \
    }                                                 \
}

#define TEST_C2_WITH_INVALID_ARGS(nG,nL,iB,str)       \
{                                                     \
    bool caught_expected_error = false;               \
    try {                                             \
      Tpetra::Map<OrdinalType> bad_map(nG,nL,iB,platform);      \
    }                                                 \
    catch (std::invalid_argument &ia) {               \
      caught_expected_error = true;                   \
      if(verbose){                                    \
        cout << str << " threw expected error:" << endl; \
      }                                               \
      if(debug){                                      \
        *fos << ia.what() << endl;                    \
      }                                               \
    }                                                 \
    if (!caught_expected_error){                      \
      ierr++;                                         \
      if(verbose){                                    \
         cout << str << " DID NOT catch expected error." << endl; \
      }                                               \
    }                                                 \
}

#define TEST_C3_WITH_INVALID_ARGS(nG,nL,eL,iB,str)    \
{                                                     \
    bool caught_expected_error = false;               \
    try {                                             \
      Tpetra::Map<OrdinalType> bad_map(nG,nL,eL,iB,platform);      \
    }                                                 \
    catch (std::invalid_argument &ia) {               \
      caught_expected_error = true;                   \
      if(verbose){                                    \
        cout << str << " threw expected error:" << endl; \
      }                                               \
      if(debug){                                      \
        *fos << ia.what() << endl;                    \
      }                                               \
    }                                                 \
    if (!caught_expected_error){                      \
      ierr++;                                         \
      if(verbose){                                    \
         cout << str << " DID NOT catch expected error." << endl; \
      }                                               \
    }                                                 \
}

#define TEST_IS_COMPATIBLE(ng1,nl1,ng2,nl2,isCompat,str)         \
{                                                      \
    Tpetra::Map<OrdinalType> m1(ng1,nl1,0,platform),   \
                             m2(ng2,nl2,0,platform);   \
    bool compat_error = false;                         \
    compat_error |= !m1.isCompatible(m1);              \
    compat_error |= !m2.isCompatible(m2);              \
    compat_error |= (m1.isCompatible(m2) != isCompat); \
    compat_error |= (m2.isCompatible(m1) != isCompat); \
    if (compat_error==false) {                         \
      if(verbose) {                                    \
         cout << str << " passed test." << endl;       \
      }                                                \
    }                                                  \
    else {                                             \
      ierr++;                                          \
      if(verbose) {                                    \
        cout << str << " FAILED test." << endl;        \
      }                                                \
    }                                                  \
}

#define TEST_IS_SAME_AS(isSameAs,str)                  \
{                                                      \
    bool issame_error = false;                         \
    /*issame_error |= !m1.isSameAs(m1);                  \
    issame_error |= !m2.isSameAs(m2);                  \
    issame_error |= (m1.isSameAs(m2) != isSameAs);     \
    issame_error |= (m2.isSameAs(m1) != isSameAs);     \
    */if (issame_error==false) {                         \
      if(verbose) {                                    \
         cout << str << " passed test." << endl;       \
      }                                                \
    }                                                  \
    else {                                             \
      ierr++;                                          \
      if(verbose) {                                    \
        cout << str << " FAILED test." << endl;        \
      }                                                \
    }                                                  \
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
  debug = debug && (myImageID == 0);

  // start the testing
  if(verbose) outputStartMessage("Map");
  int ierr = 0;

  // call the actual test routines
  ierr += unitTests<char>(verbose, debug, myImageID, numImages);
  ierr += unitTests<short int>(verbose, debug, myImageID, numImages);
  ierr += unitTests<int>(verbose, debug, myImageID, numImages);
  ierr += unitTests<long int>(verbose, debug, myImageID, numImages);

  // finish up
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  if(verbose) outputEndMessage("Map", (ierr == 0));
  return(ierr);
}


//======================================================================
template <typename OrdinalType>
int unitTests(bool verbose, bool debug, int myImageID, int numImages) {
  std::string className = "Map<" + Teuchos::OrdinalTraits<OrdinalType>::name() + ">";
  if(verbose) outputHeading("Stating unit tests for " + className);

  (void)debug;
  (void)myImageID;
  (void)numImages;

  // ======================================================================
  // code coverage section - just call functions, no testing
  // ======================================================================
  
  // create Platform and Comm
#ifdef HAVE_MPI
  Tpetra::MpiPlatform<OrdinalType>    platform(rcp(new OpaqueWrapper<MPI_Comm>(MPI_COMM_WORLD)) );
#else
  Tpetra::SerialPlatform<OrdinalType> platform;
#endif
  RCP< Comm<OrdinalType> > comm = platform.createComm();

  int ierr = 0;
  RCP<FancyOStream> fos = fancyOStream( rcp(&cout,false), std::string("  - "),1);
  // test error catching
  // the purpose of these tests is to ensure not only that errors are caught, but also that in the case that the
  // user does something stupendously incorrect (like setting a different indexBase on different images),
  // that we catch the error gracefully (all procs throw a single exception, no divergent paths, no MPI hangs)
  // that is the reason for all of the "inconsistent" checks below
  {
    // Constructor 1(numGlobal,indexBase,platform)
    // numGlobalEntries must be >= 0 and consistent across all procs
    // indexBase must be consistent across all procs
    if(verbose) cout << "Invalid arguments to constructor (numGlobal,indexBase,platform)..." << endl;
    TEST_C1_WITH_INVALID_ARGS(-1,0,"- Invalid numGlobalEntries...");
    if (numImages > 1) {
      TEST_C1_WITH_INVALID_ARGS((myImageID == 0 ? -1 : 0),0,"- Inconsistent numGlobalEntries (invalid and valid)...");
      TEST_C1_WITH_INVALID_ARGS((myImageID == 0 ?  1 : 0),0,"- Inconsistent numGlobalEntries (all valid)...");
      TEST_C1_WITH_INVALID_ARGS(0,(myImageID == 0 ? 0 : 1), "- Inconsistent indexBase...");
    }
  }
  {
    // Constructor 2(numGlobal,numLocal,indexBase,platform)
    // numGlobalEntries must be >= 0 and equal to sum numLocalEntries, or -1
    // numLocalEntries must be in [0,numGlobal], the sum must be == numGlobalEntries
    // indexBase must be consistent across all procs
    if(verbose) cout << "Invalid arguments to constructor (numGlobal,numLocal,indexBase,platform)..." << endl;
    TEST_C2_WITH_INVALID_ARGS(-2,0,0,"- Invalid numGlobalEntries...");
    TEST_C2_WITH_INVALID_ARGS(1,0,0, "- Incorrect numGlobalEntries...");
    TEST_C2_WITH_INVALID_ARGS(0,-1,0,"- Invalid numLocalEntries...");
    if (numImages > 1) {
      TEST_C2_WITH_INVALID_ARGS(numImages-2,(myImageID == 0 ? -1 : 1),0,"- Varied numLocalEntries (invalid and valid, with sum to numGlobal)");
      TEST_C2_WITH_INVALID_ARGS((myImageID == 0 ? -2 : -1),0,0,"- Inconsistent numGlobalEntries (invalid and default)...");
      TEST_C2_WITH_INVALID_ARGS((myImageID == 0 ? -2 :  0),0,0,"- Inconsistent numGlobalEntries (invalid and correct)...");
      TEST_C2_WITH_INVALID_ARGS((myImageID == 0 ? -2 :  1),0,0,"- Inconsistent numGlobalEntries (invalid and incorrect)...");
      TEST_C2_WITH_INVALID_ARGS((myImageID == 0 ? -1 :  1),0,0,"- Inconsistent numGlobalEntries (default and incorrect)...");
      TEST_C2_WITH_INVALID_ARGS((myImageID == 0 ?  1 :  0),0,0,"- Inconsistent numGlobalEntries (correct and incorrect)...");
      TEST_C2_WITH_INVALID_ARGS(0,0,(myImageID == 0 ? 0 : 1),  "- Inconsistent indexBase...");
    }
  }
  {
    // Constructor 3(numGlobal,numLocal,localEntries,indexBase,platform)
    // numGlobalEntries must be >= 0 and equal to sum numLocalEntries, or -1
    // entryList entries must be >= indexBase ??????
    // numLocalEntries must be <= entryList.size()
    // indexBase must be consistent across all procs
    if(verbose) cout << "Invalid arguments to constructor (numGlobal,numLocal,entryList,indexBase,platform)..." << endl;
    std::vector<OrdinalType> entryList(1), entryListBad(1), entryListEmpty(0);
    OrdinalType indexBase = Teuchos::OrdinalTraits<OrdinalType>::one();
    entryList[0] = myImageID + indexBase;
    entryListBad[0] = -myImageID; // \in (-numImages,0] < indexBase == 1
    TEST_C3_WITH_INVALID_ARGS(-2, 0,entryList     ,indexBase,"- Invalid numGlobalEntries...");
    TEST_C3_WITH_INVALID_ARGS( 1, 0,entryList     ,indexBase,"- Incorrect numGlobalEntries...");
    TEST_C3_WITH_INVALID_ARGS(numImages,-1,entryListEmpty,indexBase,"- Invalid numLocalEntries (negative)...");
    TEST_C3_WITH_INVALID_ARGS(numImages, 1,entryListEmpty,indexBase,"- Invalid numLocalEntries (out of range)...");
    TEST_C3_WITH_INVALID_ARGS(numImages, 1,entryListBad  ,indexBase,"- Invalid GID (less than indexBase)...");
    if (numImages > 1) {
      TEST_C3_WITH_INVALID_ARGS(numImages-2,(myImageID == 0 ? -1 : 1),entryList,indexBase,"- Varied numLocalEntries (invalid and valid, with sum to numGlobal)");
      TEST_C3_WITH_INVALID_ARGS(numImages,1,entryListBad,indexBase,"- Varied entryList (invalid and valid)");
      TEST_C3_WITH_INVALID_ARGS((myImageID == 0 ? -2 : -1),0,entryList,indexBase,"- Inconsistent numGlobalEntries (invalid and default)...");
      TEST_C3_WITH_INVALID_ARGS((myImageID == 0 ? -2 :  0),0,entryList,indexBase,"- Inconsistent numGlobalEntries (invalid and correct)...");
      TEST_C3_WITH_INVALID_ARGS((myImageID == 0 ? -2 :  1),0,entryList,indexBase,"- Inconsistent numGlobalEntries (invalid and incorrect)...");
      TEST_C3_WITH_INVALID_ARGS((myImageID == 0 ? -1 :  1),0,entryList,indexBase,"- Inconsistent numGlobalEntries (default and incorrect)...");
      TEST_C3_WITH_INVALID_ARGS((myImageID == 0 ?  1 :  0),0,entryList,indexBase,"- Inconsistent numGlobalEntries (correct and incorrect)...");
      TEST_C3_WITH_INVALID_ARGS(0,0,entryList,(myImageID == 0 ? 0 : 1),  "- Inconsistent indexBase...");
    }
  }

  // test isCompatible()
  {
    // m1.isCompatible(m2) should be true if m1 and m2 have the same number of global entries and the same number of local entries on
    // corresponding nodes
    // test the following scenarios:
    // * same number of global and local entries on all nodes
    // * same number of global entries, but different number of local entries on every node
    // * same number of global entries, but different number of local entries on some nodes
    // * different number of global entries, different number of local entries
    // 
    // for each, also:
    // test symmetry   : m1.isCompatible(m2) <=> m2.isCompatible(m1)
    // test reflexivity: m1.isCompatible(m1), m2.isCompatible(m2)
    if(verbose) cout << "Tests of Map::isCompatible()..." << endl;
    TEST_IS_COMPATIBLE(-1,myImageID,-1,myImageID,true,"- Same num global and local entries");
    TEST_IS_COMPATIBLE(-1,myImageID+1,-1,myImageID,false,"- Diff num global and local entries");
    if (numImages > 1) {
      // want different num local on every proc; map1:numLocal==[0,...,numImages-1], map2:numLocal==[1,...,numImages-1,0]
      TEST_IS_COMPATIBLE(-1,myImageID,-1,(myImageID+1)%numImages,false,"- Same num global entries, different num local entries (consistent)");
      if (numImages > 2) {
        // want different num local on a subset of procs
        // image 0 and numImages-1 get map1:numLocal==[0,numImages-1] and map2:numLocal==[numImages-1,0], the others get numLocal==myImageID
        OrdinalType mynl1, mynl2;
        if (myImageID == 0) {
          mynl1 = 0; 
          mynl2 = numImages-1;
        }
        else if (myImageID == numImages-1) {
          mynl1 = numImages-1;
          mynl2 = 0;
        }
        else {
          mynl1 = mynl2 = myImageID;
        }
        TEST_IS_COMPATIBLE(-1,mynl1,-1,mynl2,false,"- Same num global entries, different num local entries (inconsistent)");
      }
    }
  }

  // test isSameAs()
  {
    // m1.isSameAs(m2) should be true if m1 and m2 have been set equal
    // m1.isSameAs(m1) should always be true
    // m1.isSameAs(m2) is equivalent to m2.isSameAs(m1)
    // 
    // test the following scenarios:
    // * FINISH (test all divergent paths)
    // 
    if(verbose) cout << "Tests of Map::isSameAs()..." << endl;
    { 
      Tpetra::Map<OrdinalType> m1(-1,0,0,platform), m2(-1,0,0,platform);
      TEST_IS_SAME_AS(true,"- Same map (identical, empty)"); 
    }
    { 
      Tpetra::Map<OrdinalType> m1(-1,myImageID,0,platform), m2(-1,myImageID,0,platform);
      TEST_IS_SAME_AS(true,"- Same map (identical, non-empty)"); 
    }
    {
      Tpetra::Map<OrdinalType> m1(-1,myImageID,0,platform), m2(-1,myImageID+1,0,platform);
      TEST_IS_SAME_AS(false,"- Different maps (different local on each node)"); 
    }
    if (numImages > 1) {
      // FINISH: test all multi-node scenarios, esp. divergent paths
      {
        Tpetra::Map<OrdinalType> m1(-1,myImageID,0,platform), m2(-1,myImageID+(myImageID==1?1:0),0,platform);
        TEST_IS_SAME_AS(false,"- Different maps (different local on one node)"); 
      }
    }
  }

  if(verbose) {
    if(ierr == 0) {
      outputHeading("Unit tests for " + className + " passed.");
    }
    else {
      outputHeading("Unit tests for " + className + " failed: " + Tpetra::toString(ierr) + " failures.");
    }
  }
  return(ierr);
}
