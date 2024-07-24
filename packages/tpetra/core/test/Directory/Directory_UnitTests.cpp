// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Core.hpp>
#include <Tpetra_Directory.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

#include "TpetraCore_ETIHelperMacros.h"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::Map;
  using Tpetra::Directory;
  using Tpetra::LookupStatus;
  using Tpetra::IDNotPresent;
  using Tpetra::AllIDsPresent;
  using Teuchos::Array;
  using Teuchos::tuple;
  using std::sort;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Tpetra::getDefaultComm ();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, BadSize, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform map
    const GO numEntries = 2;
    RCP<M> map = rcp(new M(numEntries,0,comm));
    // create a directory
    D dir;

    Array<int> imageIDs(2);
    Array<LO> localIDs(2);
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,1)), std::invalid_argument );
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,1), localIDs(0,1)), std::invalid_argument );
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,2), localIDs(0,1)), std::invalid_argument );
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,1), localIDs(0,2)), std::invalid_argument );
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, SmallUniformContig, LO, GO )
  {
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    const LO LINV = OrdinalTraits<LO>::invalid();
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform map
    const GO numEntries = 1;
    RCP<M> map = rcp( new M(numEntries,0,comm) );
    // create a directory
    D dir;
    {
      LookupStatus stat;
      Array<int> imageIDs(numEntries);
      Array<LO>  localIDs(numEntries);
      stat = dir.getDirectoryEntries (*map, tuple<GO>(0), imageIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO>(0), imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0), localIDs );
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(numEntries+1);
      Array<LO>  localIDs(numEntries+1);
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, 1), imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, 1), imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a uniform, contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, UniformContig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    // create a uniform map
    const GO remainder = numImages/2;
    const GO numEntries = 2*numImages + remainder;
    RCP<M> map = rcp(new M(numEntries,0,comm));
    // create a directory
    D dir;
    // all GIDs
    Array<GO> allGIDs(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs[gid] = gid;
    }
    Array<int> expectedImageIDs;
    Array<LO>  expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    GO remLeft = remainder;
    for (int id = 0; id < numImages; ++id) {
      expectedImageIDs.push_back(id);
      expectedImageIDs.push_back(id);
      expectedLIDs.push_back(0);
      expectedLIDs.push_back(1);
      if (remLeft) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(2);
        --remLeft;
      }
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(numEntries);
      Array<LO> localIDs(numEntries);
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(2);
      Array<LO> localIDs(2);
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // test with a small contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, SmallContig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GO numEntries = static_cast<GO>(numImages+1);
    // the last image gets two entries, others get one
    const LO numMyEntries = (myImageID == numImages-1 ? 2 : 1);
    RCP<M> map = rcp(new M(numEntries,numMyEntries,0,comm));
    // create a directory
    D dir;
    // all GIDs
    Array<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid)
    {
      allGIDs.push_back(gid);
    }
    Array<int> expectedImageIDs;
    Array<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int id = 0; id < numImages; ++id) {
      expectedImageIDs.push_back(id);
      expectedLIDs.push_back(0);
      if (id == numImages-1) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(1);
      }
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(numEntries);
      Array<LO> localIDs(numEntries);
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(2);
      Array<LO> localIDs(2);
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // test with a contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, Contig, LO, GO )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // image i gets i+1 entries
    const LO numMyEntries = static_cast<LO>(myImageID+1);
    // number of entries is (numImages+1)*numImages/2
    const GO numEntries = static_cast<GO>((numImages*numImages+numImages)/2);
    RCP<M> map = rcp(new M(numEntries,numMyEntries,0,comm));
    // create a directory
    D dir;
    // all GIDs
    Array<GO> allGIDs(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs[gid] = gid;
    }
    Array<int> expectedImageIDs;
    Array<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int id = 0; id < numImages; ++id) {
      for (Teuchos_Ordinal num = 0; num < id+1; ++num) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(static_cast<LO>(num));
      }
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(numEntries);
      Array<LO> localIDs(numEntries);
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(2);
      Array<LO> localIDs(2);
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  // test with a non-contiguous map
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Directory, NonContig, LO, GO )
  {
    using std::endl;
    typedef Map<LO,GO> M;
    typedef Directory<LO,GO> D;

    const LO LINV = OrdinalTraits<LO>::invalid();
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // number of entries is 3*numImages
    // we will stripe the GIDs across images
    const GO numEntries = static_cast<GO>(3*numImages);

    out << "Creating Map" << endl;
    RCP<M> map = rcp (new M (numEntries, tuple<GO> (myImageID, myImageID+numImages, myImageID+2*numImages), 0, comm));

    out << "Creating Directory" << endl;
    D dir;

    out << "Create array of GIDs to look up" << endl;
    Array<GO> allGIDs;
    allGIDs.reserve(numEntries);
    for (GO gid = 0; gid < numEntries; ++gid) {
      allGIDs.push_back(gid);
    }
    Array<int> expectedImageIDs;
    Array<LO> expectedLIDs;
    // expected image IDs and LIDs
    expectedImageIDs.reserve(numEntries);
    expectedLIDs.reserve(numEntries);
    for (int i = 0; i < 3; ++i) {
      for (int id = 0; id < numImages; ++id) {
        expectedImageIDs.push_back(id);
        expectedLIDs.push_back(i);
      }
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(numEntries);
      Array<LO> localIDs(numEntries);
      stat = dir.getDirectoryEntries (*map, allGIDs,imageIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      stat = dir.getDirectoryEntries (*map, allGIDs, imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, AllIDsPresent );
      TEST_COMPARE_ARRAYS( expectedImageIDs, imageIDs );
      TEST_COMPARE_ARRAYS( expectedLIDs, localIDs );
    }
    {
      LookupStatus stat;
      Array<int> imageIDs(2);
      Array<LO> localIDs(2);
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries), imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallUniformContig, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, UniformContig, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, SmallContig, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, Contig, LO, GO )     \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, NonContig, LO, GO )  \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Directory, BadSize, LO, GO )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

} // namespace (anonymous)



