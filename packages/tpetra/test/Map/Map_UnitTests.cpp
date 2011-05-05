#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Map.hpp"

// FINISH: add testing of operator==, operator!=, operator=, copy construct
// put these into test_same_as and test_is_compatible

namespace {

  using Teuchos::Array;
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

  bool testMpi = true;
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
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignord and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  // 


  ////
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


  ////
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


  ////
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
    if (m != Teuchos::null) {
      TEST_EQUALITY( m->getMinAllGlobalIndex(), (size_t)0 );
      TEST_EQUALITY( m->getMaxAllGlobalIndex(), (size_t)numImages-1 );
    }
  }
#endif


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, compatabilityTests, LO, GO )
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
    M map(numGlobalEntries,indexBase,comm);

    TEST_EQUALITY_CONST(map.isContiguous(), true);
    TEST_EQUALITY_CONST(map.isDistributed(), numImages > 1);
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobalEntries);
    TEST_EQUALITY_CONST(map.getNodeNumElements(), 2);
    TEST_EQUALITY_CONST(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), indexBase);
    TEST_EQUALITY_CONST(map.getMaxLocalIndex(), 1);
    TEST_EQUALITY_CONST(map.getMinGlobalIndex(), myGlobal[indexBase]);
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


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor1, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor2, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor3, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, compatabilityTests, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, sameasTests, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ContigUniformMap, LO, GO )

    UNIT_TEST_GROUP_ORDINAL( char , int )
    UNIT_TEST_GROUP_ORDINAL( int , int )

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor1, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor2, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, invalidConstructor3, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, compatabilityTests, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, sameasTests, LO, GO ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ContigUniformMap, LO, GO )

    // UNIT_TEST_GROUP_ORDINAL(char , int)

    UNIT_TEST_GROUP_ORDINAL(int , int)

    // typedef short int ShortInt;
    // UNIT_TEST_GROUP_ORDINAL(ShortInt, int)

    // typedef long int LongInt;
    // UNIT_TEST_GROUP_ORDINAL(int , LongInt)

#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      // typedef long long int LongLongInt;
      // UNIT_TEST_GROUP_ORDINAL(char , LongLongInt)
      // UNIT_TEST_GROUP_ORDINAL(int , LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
