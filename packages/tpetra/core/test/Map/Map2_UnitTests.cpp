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

#include "TpetraNew_Map.hpp"
#include "Tpetra_TestingUtilities.hpp"
#include <type_traits> // std::is_same

namespace {

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Teuchos::tuple;
  using Teuchos::OrdinalTraits;
  using std::endl;

  using LO = TpetraNew::Map::local_ordinal_type;
  using GO = TpetraNew::Map::global_ordinal_type;  

#define TEST_IS_COMPATIBLE(m1,m2,is_compat)               \
{                                                         \
  Teuchos::OSTab tabCompat0 (out);                        \
  out << "Expect " << (is_compat ? "" : "NOT ") << "compatible" << std::endl; \
  Teuchos::OSTab tabCompat1 (out); \
  out << "Is m1 compatible with itself?" << std::endl;  \
  TEST_EQUALITY_CONST(m1.isCompatible(m1), true);       \
  out << "Is m2 compatible with itself?" << std::endl;  \
  TEST_EQUALITY_CONST(m2.isCompatible(m2), true);       \
  out << "Is m1 compatible with m2?" << std::endl;      \
  TEST_EQUALITY_CONST(m1.isCompatible(m2), is_compat);  \
  out << "Is m2 compatible with m1?" << std::endl;      \
  TEST_EQUALITY_CONST(m2.isCompatible(m1), is_compat);  \
}

#define TEST_IS_SAME_AS(m1,m2,is_sameas)               \
{                                                      \
  out << "Expect " << (is_sameas ? "" : "NOT ") << "same" << std::endl; \
    TEST_EQUALITY_CONST(m1.isSameAs(m1), true);        \
    TEST_EQUALITY_CONST(m2.isSameAs(m2), true);        \
    TEST_EQUALITY_CONST(m1.isSameAs(m2), is_sameas);   \
    TEST_EQUALITY_CONST(m2.isSameAs(m1), is_sameas);   \
}

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  //
  // UNIT TESTS
  //

  // This test exercises unsigned GlobalOrdinal.  We don't intend to
  // support unsigned GlobalOrdinal going forward.  However, the test
  // should work regardless of the GlobalOrdinal type.
  TEUCHOS_UNIT_TEST( Map2, RogersUnsignedGOBugVerification )
  {
    using map_type = TpetraNew::Map;
    
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 2) return;
    const int myImageID = comm->getRank();
    const GO INVALID = OrdinalTraits<GO>::invalid();
    RCP<map_type> m;
    TEST_NOTHROW( m = rcp(new map_type(INVALID, tuple<GO>(myImageID), 0, comm)) );
    TEST_ASSERT( m.get () != nullptr );
    if (m.get () != nullptr) {
      TEST_EQUALITY( m->getMinAllGlobalIndex(),
		     static_cast<GO> (0) );
      TEST_EQUALITY( m->getMaxAllGlobalIndex(),
		     static_cast<GO> (numImages-1) );
    }

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST( Map2, compatibilityTests )
  {
    using std::endl;
    typedef TpetraNew::Map M;

    out << "Test: Map, compatibilityTests" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GO INVALID = OrdinalTraits<GO>::invalid();

    // test isCompatible()
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
    {
      out << "Contiguous nonuniform ctor, same local GID count globally" << endl;
      M m1(INVALID,myImageID,0,comm),
        m2(INVALID,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, true );
    }
    {
      out << "Contiguous nonuniform ctor, different local and global GID counts" << endl;
      M m1(INVALID,myImageID+1,0,comm),
        m2(INVALID,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, false);
    }
    if (numImages > 1) {
      // want different num local on every proc; map1:numLocal==[0,...,numImages-1], map2:numLocal==[1,...,numImages-1,0]
      {
        out << "Contiguous nonuniform ctor, same global GID count, "
            << "different local GID counts" << endl;
        M m1(INVALID,myImageID,0,comm),
          m2(INVALID,(myImageID+1)%numImages,0,comm);
        out << "myImageID = " << myImageID
            << ", (myImageID+1) % numImages = "
            << ((myImageID+1) % numImages) << endl;
        Teuchos::OSTab tab1 (out);
        out << "m1.getGlobalNumElements() = " << m1.getGlobalNumElements () << endl
            << "m2.getGlobalNumElements() = " << m2.getGlobalNumElements () << endl
            << "m1.getNodeNumElements() = " << m1.getNodeNumElements () << endl
            << "m2.getNodeNumElements() = " << m2.getNodeNumElements () << endl;
        TEST_IS_COMPATIBLE( m1, m2, false);
      }
      if (numImages > 2) {
        // want different num local on a subset of procs
        // image 0 and numImages-1 get map1:numLocal==[0,numImages-1] and map2:numLocal==[numImages-1,0], the others get numLocal==myImageID
        LO mynl1, mynl2;
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
        {
          out << "Contiguous nonuniform ctor, same global GID count, "
              << "different local GID counts on subset of processes" << endl;
          M m1(INVALID,mynl1,0,comm),
            m2(INVALID,mynl2,0,comm);
          TEST_IS_COMPATIBLE( m1, m2, false);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST( Map2, sameAsTests )
  {
    typedef TpetraNew::Map M;

    out << "Test: Map, sameasTests" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GO INVALID = OrdinalTraits<GO>::invalid();
    {
      M m1(INVALID,0,0,comm),
        m2(INVALID,0,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(INVALID,myImageID,0,comm),
        m2(INVALID,myImageID,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(INVALID,myImageID,0,comm),
        m2(INVALID,myImageID+1,0,comm);
      TEST_IS_SAME_AS(m1, m2, false);
    }
    if (numImages > 1) {
      // FINISH: test all multi-node scenarios, esp. divergent paths
      {
        M m1(INVALID,myImageID,0,comm),
          m2(INVALID,myImageID+(myImageID==1?1:0),0,comm);
        TEST_IS_SAME_AS(m1, m2, false);
      }
    }
  }

  TEUCHOS_UNIT_TEST( Map2, ContigUniformMap )
  {
    typedef TpetraNew::Map M;

    out << "Test: Map, ContigUniformMap" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with two entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*2, myImageID*2+1) );
    Array<LO>  myLocal( tuple<LO>(0,1) );

    const GO numGlobalEntries = GO (numImages) * GO (2);
    const GO indexBase = 0;
    const LO localIndexBase = 0;
    M map(numGlobalEntries,indexBase,comm);

    TEST_EQUALITY_CONST(map.isContiguous(), true);
    TEST_EQUALITY_CONST(map.isDistributed(), numImages > 1);
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobalEntries);
    TEST_EQUALITY_CONST(map.getNodeNumElements(), LO (2));
    TEST_EQUALITY_CONST(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), localIndexBase);
    TEST_EQUALITY_CONST(map.getMaxLocalIndex(), LO (1));
    TEST_EQUALITY_CONST(map.getMinGlobalIndex(), myGlobal[0]);
    TEST_EQUALITY_CONST(map.getMaxGlobalIndex(), myGlobal[1]);
    TEST_EQUALITY_CONST(map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST(map.getMaxAllGlobalIndex(),
			static_cast<GO>(numGlobalEntries-1));
    TEST_EQUALITY( map.getLocalElement(myGlobal[0]), myLocal[0] );
    TEST_EQUALITY( map.getLocalElement(myGlobal[1]), myLocal[1] );
    TEST_EQUALITY( map.getGlobalElement(myLocal[0]), myGlobal[0] );
    TEST_EQUALITY( map.getGlobalElement(myLocal[1]), myGlobal[1] );
    TEST_EQUALITY( map.getLocalElement(numGlobalEntries), OrdinalTraits<LO>::invalid() );
    TEST_EQUALITY( map.getGlobalElement(2),               OrdinalTraits<GO>::invalid() );
    TEST_EQUALITY( map.getLocalElement(numGlobalEntries-1), myImageID == numImages-1 ? 1 : OrdinalTraits<LO>::invalid() );
    TEST_COMPARE_ARRAYS( map.getNodeElementList(), myGlobal);
    TEST_EQUALITY_CONST( map.isNodeLocalElement(0), true );
    TEST_EQUALITY_CONST( map.isNodeLocalElement(1), true );
    TEST_EQUALITY_CONST( map.isNodeLocalElement(2), false ); // just try a couple
    TEST_EQUALITY_CONST( map.isNodeLocalElement(3), false );
    for (GO i=0; i < static_cast<GO>(numGlobalEntries); ++i) {
      if (std::find(myGlobal.begin(),myGlobal.end(),i) == myGlobal.end()) {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), false );
      }
      else {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), true );
      }
    }

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

TEUCHOS_UNIT_TEST( Map2, nonTrivialIndexBase )
  {
    typedef TpetraNew::Map Map;

    out << "Test: Map, nonTrivialIndexBase" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const LO numLocal  = 5;
    const GO numGlobal = GO (numImages) * GO (numLocal);
    const GO indexBase = 10;

    Map map (numGlobal, indexBase, comm);
    //
    Array<GO> expectedGIDs(numLocal);
    expectedGIDs[0] = indexBase + myImageID*numLocal;
    for (LO i=1; i<numLocal; ++i) {
      expectedGIDs[i] = expectedGIDs[i-1]+1;
    }
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getNodeNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), 0);
    TEST_EQUALITY(map.getMaxLocalIndex(), numLocal-1);
    TEST_EQUALITY(map.getMinGlobalIndex(), expectedGIDs[0]);
    TEST_EQUALITY(map.getMaxGlobalIndex(),
		  static_cast<GO>(expectedGIDs[0]+numLocal-1) );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY(map.getGlobalElement(0), expectedGIDs[0]);
    TEST_EQUALITY_CONST(static_cast<GO> (map.getMaxAllGlobalIndex()),
			static_cast<GO> (indexBase+numGlobal-1));
    ArrayView<const GO> glist = map.getNodeElementList();
    TEST_COMPARE_ARRAYS( map.getNodeElementList(), expectedGIDs);

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

TEUCHOS_UNIT_TEST( Map2, indexBaseAndAllMin )
  {
    typedef TpetraNew::Map Map;

    out << "Test: Map, indexBaseAndAllMin" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const LO numLocal = 5;
    const GO numGlobal = GO (numImages) * GO (numLocal);
    const GO indexBase = 0;
    const GO actualBase = 1;
    //
    Array<GO> GIDs(numLocal);
    GIDs[0] = actualBase + myImageID*numLocal;
    for (LO i=1; i<numLocal; ++i) {
      GIDs[i] = GIDs[i-1]+1;
    }

    Map map (numGlobal, GIDs (), indexBase, comm);
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getNodeNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), 0);
    TEST_EQUALITY(map.getMaxLocalIndex(), numLocal-1);
    TEST_EQUALITY(map.getMinGlobalIndex(), GIDs[0]);
    TEST_EQUALITY(map.getMaxGlobalIndex(), static_cast<GO>(GIDs[0]+numLocal-1) );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), actualBase);
    TEST_EQUALITY(map.getGlobalElement(0), GIDs[0]);
    TEST_EQUALITY_CONST(static_cast<GO> (map.getMaxAllGlobalIndex ()),
			static_cast<GO> (actualBase+numGlobal-1));
    ArrayView<const GO> glist = map.getNodeElementList();
    TEST_COMPARE_ARRAYS( map.getNodeElementList(), GIDs);

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

TEUCHOS_UNIT_TEST( Map2, ZeroLocalElements )
  {
    typedef TpetraNew::Map M;

    out << "Test: Map, ZeroLocalElements" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int rank = comm->getRank();

    // Create maps with zero elements on all but the first processor
    Array<GO>  elem_list;
    if (rank == 0)
      elem_list.push_back(0);
    M contig_uniform(1, 0, comm);
    M contig_non_uniform(1, elem_list.size(), 0, comm);
    M non_contig(1, elem_list, 0, comm);

    // Check LID
    LO lid_expected = rank == 0 ? LO (0) : OrdinalTraits<LO>::invalid();
    TEST_EQUALITY( contig_uniform.getLocalElement(0), lid_expected );
    TEST_EQUALITY( contig_non_uniform.getLocalElement(0), lid_expected );
    TEST_EQUALITY( non_contig.getLocalElement(0), lid_expected );

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST( Map, AllDefaultTemplateParameters )
  {
    out << "Test: Map, AllDefaultTemplateParameters" << std::endl;
    Teuchos::OSTab tab0 (out);
    static_assert (sizeof (GO) >= sizeof (LO),
		   "Map::global_ordinal_type must have size no "
		   "less than that of Map::local_ordinal_type.");
  }

} // namespace (anonymous)


