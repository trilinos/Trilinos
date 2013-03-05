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

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_Map.hpp>

// FINISH: add testing of operator==, operator!=, operator=, copy construct
// put these into test_same_as and test_is_compatible

namespace {

  using Teuchos::null;

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using Tpetra::createUniformContigMapWithNode;

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using Tpetra::Map;
  using Tpetra::global_size_t;
  using Tpetra::DefaultPlatform;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  double errorTolSlack = 1e+1;

#define TEST_IS_COMPATIBLE(m1,m2,is_compat)               \
{                                                         \
    TEST_EQUALITY_CONST(m1.isCompatible(m1), true);       \
    TEST_EQUALITY_CONST(m2.isCompatible(m2), true);       \
    TEST_EQUALITY_CONST(m1.isCompatible(m2), is_compat);  \
    TEST_EQUALITY_CONST(m2.isCompatible(m1), is_compat);  \
}

#define TEST_IS_SAME_AS(m1,m2,is_sameas)               \
{                                                      \
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
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //

#if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) && defined(HAVE_TPETRA_ENABLE_SS_TESTING) && defined(HAVE_TPETRA_MPI)
  ////
  TEUCHOS_UNIT_TEST( Map, RogersUnsignedGOBugVerification )
  {
    typedef Map<int,size_t> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 2) return;
    const int myImageID = comm->getRank();
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    RCP<M> m;
    TEST_NOTHROW( m = rcp(new M(GSTI, tuple<size_t>(myImageID), 0, comm)) );
    if (m != null) {
      TEST_EQUALITY( m->getMinAllGlobalIndex(), (size_t)0 );
      TEST_EQUALITY( m->getMaxAllGlobalIndex(), (size_t)numImages-1 );
    }
  }
#endif

  // This test may only pass in a debug build (HAVE_TPETRA_DEBUG).
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, invalidConstructor1, LO, GO )
  {
    typedef Map<LO,GO> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    // bad constructor calls: (num global elements, index base)
    TEST_THROW(M map(GSTI,0,comm), std::invalid_argument);
    if (numImages > 1) {
      TEST_THROW(M map((myImageID == 0 ? GSTI : 0),0,comm), std::invalid_argument);
      TEST_THROW(M map((myImageID == 0 ?  1 : 0),0,comm), std::invalid_argument);
      TEST_THROW(M map(0,(myImageID == 0 ? 0 : 1), comm), std::invalid_argument);
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // This test may only pass in a debug build (HAVE_TPETRA_DEBUG).
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, invalidConstructor2, LO, GO )
  {
    typedef Map<LO,GO> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    // bad constructor calls: (num global elements, num local elements, index base)
    TEST_THROW(M map(1,0,0, comm),  std::invalid_argument);
    if (numImages > 1) {
      TEST_THROW(M map((myImageID == 0 ? GSTI :  1),0,0,comm), std::invalid_argument);
      TEST_THROW(M map((myImageID == 0 ?  1 :  0),0,0,comm), std::invalid_argument);
      TEST_THROW(M map(0,0,(myImageID == 0 ? 0 : 1),comm), std::invalid_argument);
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  // This test may only pass in a debug build (HAVE_TPETRA_DEBUG).
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, invalidConstructor3, LO, GO )
  {
    typedef Map<LO,GO> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    // bad constructor calls: (num global, entry list, index base)
    TEST_THROW(M map(numImages, tuple<GO>(-myImageID), 1, comm), std::invalid_argument); // GID less than iB
    if (numImages > 1) {
      TEST_THROW(M map( 1, tuple<GO>(myImageID+1), 1, comm), std::invalid_argument);    // nG != sum nL
      TEST_THROW(M map((myImageID == 0 ? GSTI :  0),tuple<GO>(myImageID+1),1, comm), std::invalid_argument);
      TEST_THROW(M map(0, tuple<GO>(myImageID+1), (myImageID == 0 ? 0 : 1), comm), std::invalid_argument);
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, compatibilityTests, LO, GO )
  {
    typedef Map<LO,GO> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
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
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, true );
    }
    {
      M m1(GSTI,myImageID+1,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, false);
    }
    if (numImages > 1) {
      // want different num local on every proc; map1:numLocal==[0,...,numImages-1], map2:numLocal==[1,...,numImages-1,0]
      {
        M m1(GSTI,myImageID,0,comm),
          m2(GSTI,(myImageID+1)%numImages,0,comm);
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
          M m1(GSTI,mynl1,0,comm),
            m2(GSTI,mynl2,0,comm);
          TEST_IS_COMPATIBLE( m1, m2, false);
        }
      }
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, sameasTests, LO, GO )
  {
    typedef Map<LO,GO> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    {
      M m1(GSTI,0,0,comm),
        m2(GSTI,0,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID+1,0,comm);
      TEST_IS_SAME_AS(m1, m2, false);
    }
    if (numImages > 1) {
      // FINISH: test all multi-node scenarios, esp. divergent paths
      {
        M m1(GSTI,myImageID,0,comm),
          m2(GSTI,myImageID+(myImageID==1?1:0),0,comm);
        TEST_IS_SAME_AS(m1, m2, false);
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, ContigUniformMap, LO, GO )
  {
    typedef Map<LO,GO> M;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with two entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*2, myImageID*2+1) );
    Array<LO>  myLocal( tuple<LO>(0,1) );

    const size_t numGlobalEntries = numImages*2;
    const GO indexBase = 0;
    const LO localIndexBase = 0;
    M map(numGlobalEntries,indexBase,comm);

    TEST_EQUALITY_CONST(map.isContiguous(), true);
    TEST_EQUALITY_CONST(map.isDistributed(), numImages > 1);
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobalEntries);
    TEST_EQUALITY_CONST(map.getNodeNumElements(), 2);
    TEST_EQUALITY_CONST(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), localIndexBase);
    TEST_EQUALITY_CONST(map.getMaxLocalIndex(), 1);
    TEST_EQUALITY_CONST(map.getMinGlobalIndex(), myGlobal[0]);
    TEST_EQUALITY_CONST(map.getMaxGlobalIndex(), myGlobal[1]);
    TEST_EQUALITY_CONST(map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST(map.getMaxAllGlobalIndex(), as<GO>(numGlobalEntries-1));
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
    for (GO i=0; i < as<GO>(numGlobalEntries); ++i) {
      if (find(myGlobal.begin(),myGlobal.end(),i) == myGlobal.end()) {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), false );
      }
      else {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), true );
      }
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, nonTrivialIndexBase, LO, GO )
  {
    typedef Map<LO,GO> Map;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t        numLocal  = 5;
    const global_size_t numGlobal = numImages*numLocal;
    const GO indexBase = 10;
    RCP<Map> map = rcp(new Map(numGlobal,indexBase,comm));
    //
    Array<GO> expectedGIDs(numLocal);
    expectedGIDs[0] = indexBase + myImageID*numLocal;
    for (size_t i=1; i<numLocal; ++i) {
      expectedGIDs[i] = expectedGIDs[i-1]+1;
    }
    //
    TEST_EQUALITY(map->getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map->getNodeNumElements(), numLocal);
    TEST_EQUALITY(map->getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map->getMinLocalIndex(), 0);
    TEST_EQUALITY(map->getMaxLocalIndex(), numLocal-1);
    TEST_EQUALITY(map->getMinGlobalIndex(), expectedGIDs[0]);
    TEST_EQUALITY(map->getMaxGlobalIndex(), as<GO>(expectedGIDs[0]+numLocal-1) );
    TEST_EQUALITY(map->getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY(map->getGlobalElement(0), expectedGIDs[0]);
    TEST_EQUALITY_CONST((global_size_t)map->getMaxAllGlobalIndex(), indexBase+numGlobal-1);
    ArrayView<const GO> glist = map->getNodeElementList();
    TEST_COMPARE_ARRAYS( map->getNodeElementList(), expectedGIDs);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, indexBaseAndAllMin, LO, GO )
  {
    typedef Map<LO,GO> Map;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t        numLocal  = 5;
    const global_size_t numGlobal = numImages*numLocal;
    const GO indexBase = 0;
    const GO actualBase = 1;
    //
    Array<GO> GIDs(numLocal);
    GIDs[0] = actualBase + myImageID*numLocal;
    for (size_t i=1; i<numLocal; ++i) {
      GIDs[i] = GIDs[i-1]+1;
    }
    RCP<Map> map = rcp(new Map(numGlobal,GIDs(),indexBase,comm));
    //
    TEST_EQUALITY(map->getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map->getNodeNumElements(), numLocal);
    TEST_EQUALITY(map->getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map->getMinLocalIndex(), 0);
    TEST_EQUALITY(map->getMaxLocalIndex(), numLocal-1);
    TEST_EQUALITY(map->getMinGlobalIndex(), GIDs[0]);
    TEST_EQUALITY(map->getMaxGlobalIndex(), as<GO>(GIDs[0]+numLocal-1) );
    TEST_EQUALITY(map->getMinAllGlobalIndex(), actualBase);
    TEST_EQUALITY(map->getGlobalElement(0), GIDs[0]);
    TEST_EQUALITY_CONST((global_size_t)map->getMaxAllGlobalIndex(), actualBase+numGlobal-1);
    ArrayView<const GO> glist = map->getNodeElementList();
    TEST_COMPARE_ARRAYS( map->getNodeElementList(), GIDs);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Map, NodeConversion, LO, GO, N2 )
  {
    typedef typename Kokkos::DefaultNode::DefaultNodeType N1;
    typedef Map<LO,GO,N1> Map1;
    typedef Map<LO,GO,N2> Map2;
    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();
    const size_t        numLocal  = 10;
    const global_size_t numGlobal = numImages*numLocal;

    RCP<N1> n1 = getNode<N1>();
    RCP<N2> n2 = getNode<N2>();

    // create a contiguous uniform distributed map with numLocal entries per node
    RCP<const Map1> map1 = createUniformContigMapWithNode<LO,GO>(numGlobal,comm,n1);
    RCP<const Map2> map2 = map1->clone(n2);
    TEST_EQUALITY( map2->getNode(), n2 );
    RCP<const Map1> map1b = map2->clone(n1);
    TEST_EQUALITY( map1b->getNode(), n1 );
    TEST_EQUALITY_CONST( map1->isCompatible(*map1b), true );
    TEST_EQUALITY_CONST( map1->isSameAs(*map1b), true );
  }


  //
  // INSTANTIATIONS
  //

#ifdef HAVE_TPETRA_DEBUG
  // all ordinals, default node
#  define UNIT_TEST_GROUP( LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor1, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor2, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor3, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, compatibilityTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, sameasTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, nonTrivialIndexBase, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, indexBaseAndAllMin, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ContigUniformMap, LO, GO )
#else
  // all ordinals, default node
#  define UNIT_TEST_GROUP( LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, compatibilityTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, sameasTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, nonTrivialIndexBase, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, indexBaseAndAllMin, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ContigUniformMap, LO, GO )
#endif // HAVE_TPETRA_DEBUG

#define NC_TESTS(N) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, NodeConversion, int, int, N )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

  TPETRA_INSTANTIATE_N(NC_TESTS)

}
