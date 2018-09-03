/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_Core.hpp>
#include <TpetraNew_Directory.hpp>
#include <TpetraNew_Map.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>
#include <numeric>

namespace {

  using Tpetra::LookupStatus;
  using Tpetra::IDNotPresent;
  using Tpetra::AllIDsPresent;
  using Teuchos::Array;
  using Teuchos::broadcast;
  using Teuchos::Comm;  
  using Teuchos::OrdinalTraits;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::endl;
  using std::sort;
  
  using M = TpetraNew::Map;
  using D = TpetraNew::Directory;
  using LO = M::local_ordinal_type;
  using GO = M::global_ordinal_type;
  
  bool testMpi = true;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Tpetra::getDefaultComm ();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST( Directory_new, BadSize )
  {
    auto comm = getDefaultComm ();
    const GO numEntries = 2;
    auto map = rcp (new M (numEntries, 0, comm));
    D dir;

    Array<int> imageIDs(2);
    Array<LO> localIDs(2);
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,1)), std::invalid_argument );
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,1), localIDs(0,1)), std::invalid_argument );
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,2), localIDs(0,1)), std::invalid_argument );
    TEST_THROW( dir.getDirectoryEntries(*map, tuple<GO>(0,1), imageIDs(0,1), localIDs(0,2)), std::invalid_argument );
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a uniform, contiguous map of constant size
  TEUCHOS_UNIT_TEST( Directory_new, SmallUniformContig )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    auto comm = getDefaultComm ();
    const GO numEntries = 1;
    auto map = rcp (new M (numEntries, 0, comm));
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
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a uniform, contiguous map
  TEUCHOS_UNIT_TEST( Directory, UniformContig )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    auto comm = getDefaultComm();
    const int numImages = comm->getSize();
    const GO remainder = numImages/2;
    const GO numEntries = 2*numImages + remainder;
    auto map = rcp (new M (numEntries, 0, comm));
    D dir;
    Array<GO> allGIDs(numEntries);
    std::iota (allGIDs.begin (), allGIDs.end (), GO (0));
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
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries),
				      imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries),
				      imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM,
			success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a small contiguous map
  TEUCHOS_UNIT_TEST( Directory_new, SmallContig )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    auto comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GO numEntries = static_cast<GO>(numImages+1);
    // the last process gets two entries, others get one
    const LO numMyEntries = (myImageID == numImages-1 ? 2 : 1);
    auto map = rcp (new M (numEntries, numMyEntries, 0, comm));
    D dir;
    Array<GO> allGIDs (numEntries);
    std::iota (allGIDs.begin (), allGIDs.end (), GO (0));
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
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries),
				      imageIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      stat = dir.getDirectoryEntries (*map, tuple<GO> (0, numEntries),
				      imageIDs, localIDs);
      TEST_EQUALITY_CONST( stat, IDNotPresent );
      TEST_COMPARE_ARRAYS( tuple<int>(0,-1), imageIDs );
      TEST_COMPARE_ARRAYS( tuple<LO>(0,LINV), localIDs );
    }
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM,
			success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a contiguous map
  TEUCHOS_UNIT_TEST( Directory_new, Contig )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    auto comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // image i gets i+1 entries
    const LO numMyEntries = static_cast<LO>(myImageID+1);
    // number of entries is (numImages+1)*numImages/2
    const GO numEntries =
      (GO (numImages) * GO (numImages) + GO (numImages)) / GO (2);
    auto map = rcp (new M (numEntries, numMyEntries, 0, comm));
    D dir;
    Array<GO> allGIDs(numEntries);
    std::iota (allGIDs.begin (), allGIDs.end (), GO (0));
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
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // test with a non-contiguous map
  TEUCHOS_UNIT_TEST( Directory_new, NonContig )
  {
    const LO LINV = OrdinalTraits<LO>::invalid();
    auto comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // number of entries is 3*numImages
    // we will stripe the GIDs across images
    const GO numEntries = static_cast<GO>(3*numImages);

    out << "Creating Map" << endl;
    auto map = rcp (new M (numEntries, tuple<GO> (myImageID, myImageID+numImages, myImageID+2*numImages), 0, comm));

    out << "Creating Directory" << endl;
    D dir;

    out << "Create array of GIDs to look up" << endl;
    Array<GO> allGIDs (numEntries);
    std::iota (allGIDs.begin (), allGIDs.end (), GO (0));
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
    // All procs fail if any process fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM,
			success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

} // namespace (anonymous)



