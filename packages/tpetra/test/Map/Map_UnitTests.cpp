#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Map.hpp"

// FINISH: add testing of operator==, operator!=, operator=, copy construct
// put these into test_same_as and test_is_compatible

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Tpetra::Map;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
   }

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

  template<class Ordinal>
  RCP<const Platform<Ordinal> > getDefaultPlatform()
  {
    if (testMpi) {
      return DefaultPlatform<Ordinal>::getPlatform();
    }
    return rcp(new Tpetra::SerialPlatform<Ordinal>());
  }

  //
  // UNIT TESTS
  // 


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Map, invalidConstructor1, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // bad constructor calls: (global entries, index base)
    TEST_THROW(Map<Ordinal> map(-1,0,platform), std::invalid_argument);
    if (numImages > 1) {
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -1 : 0),0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ?  1 : 0),0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map(0,(myImageID == 0 ? 0 : 1), platform), std::invalid_argument);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Map, invalidConstructor2, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // bad constructor calls: (global entries, local entries, index base, )
    TEST_THROW(Map<Ordinal> map(-2,0,0,platform), std::invalid_argument);
    TEST_THROW(Map<Ordinal> map(1,0,0,platform),  std::invalid_argument);
    TEST_THROW(Map<Ordinal> map(0,-1,0,platform), std::invalid_argument);
    if (numImages > 1) {
      TEST_THROW(Map<Ordinal> map(numImages-2,(myImageID == 0 ? -1 : 1),0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -2 : -1),0,0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -2 :  0),0,0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -2 :  1),0,0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -1 :  1),0,0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ?  1 :  0),0,0,platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map(0,0,(myImageID == 0 ? 0 : 1),platform), std::invalid_argument);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Map, invalidConstructor3, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // bad constructor calls: (num global, entry list, index base)
    TEST_THROW(Map<Ordinal> map(-2, arcp(rcp(new vector<Ordinal>(1,myImageID+1))).getConst()     ,1, platform), std::invalid_argument);    // nG not valid
    TEST_THROW(Map<Ordinal> map(numImages, arcp(rcp(new vector<Ordinal>(1,-myImageID))).getConst()  ,1, platform), std::invalid_argument); // GID less than iB
    if (numImages > 1) {
      TEST_THROW(Map<Ordinal> map( 1, arcp(rcp(new vector<Ordinal>(1,myImageID+1))).getConst()     ,1, platform), std::invalid_argument);    // nG != sum nL
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -2 : -1),arcp(rcp(new vector<Ordinal>(1,myImageID+1))).getConst(),1, platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -2 :  0),arcp(rcp(new vector<Ordinal>(1,myImageID+1))).getConst(),1, platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map((myImageID == 0 ? -1 :  0),arcp(rcp(new vector<Ordinal>(1,myImageID+1))).getConst(),1, platform), std::invalid_argument);
      TEST_THROW(Map<Ordinal> map(0,arcp(rcp(new vector<Ordinal>(1,myImageID+1))).getConst(),(myImageID == 0 ? 0 : 1), platform), std::invalid_argument);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Map, compatabilityTests, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
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
    TEST_IS_COMPATIBLE( Map<Ordinal>(-1,myImageID,0,platform), Map<Ordinal>(-1,myImageID,0,platform), true );
    TEST_IS_COMPATIBLE( Map<Ordinal>(-1,myImageID+1,0,platform), Map<Ordinal>(-1,myImageID,0,platform), false);
    if (numImages > 1) {
      // want different num local on every proc; map1:numLocal==[0,...,numImages-1], map2:numLocal==[1,...,numImages-1,0]
      TEST_IS_COMPATIBLE( Map<Ordinal>(-1,myImageID,0,platform), Map<Ordinal>(-1,(myImageID+1)%numImages,0,platform), false);
      if (numImages > 2) {
        // want different num local on a subset of procs
        // image 0 and numImages-1 get map1:numLocal==[0,numImages-1] and map2:numLocal==[numImages-1,0], the others get numLocal==myImageID
        Ordinal mynl1, mynl2;
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
        TEST_IS_COMPATIBLE( Map<Ordinal>(-1,mynl1,0,platform), Map<Ordinal>(-1,mynl2,0,platform), false);
      }
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Map, sameasTests, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    TEST_IS_SAME_AS(Map<Ordinal>(-1,0,0,platform), Map<Ordinal>(-1,0,0,platform), true);
    TEST_IS_SAME_AS(Map<Ordinal>(-1,myImageID,0,platform), Map<Ordinal>(-1,myImageID,0,platform), true);
    TEST_IS_SAME_AS(Map<Ordinal>(-1,myImageID,0,platform), Map<Ordinal>(-1,myImageID+1,0,platform), false);
    if (numImages > 1) {
      // FINISH: test all multi-node scenarios, esp. divergent paths
      TEST_IS_SAME_AS(Map<Ordinal>(-1,myImageID,0,platform), Map<Ordinal>(-1,myImageID+(myImageID==1?1:0),0,platform), false);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Map, ContigUniformMap, Ordinal )
  {
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with two entries per node
    // this map will have the following entries:
    vector<Ordinal> myLocal, myGlobal;
    myLocal.push_back(as<Ordinal>(0));
    myLocal.push_back(as<Ordinal>(1));
    myGlobal.push_back(as<Ordinal>(myImageID*2));
    myGlobal.push_back(as<Ordinal>(myImageID*2+1));

    const Ordinal numGlobalEntries = as<Ordinal>(numImages*2);
    const Ordinal indexBase = as<Ordinal>(0);
    Map<Ordinal> map(numGlobalEntries,indexBase,platform);

    TEST_EQUALITY_CONST(map.isContiguous(), true);
    TEST_EQUALITY_CONST(map.isDistributed(), numImages > 1);
    TEST_EQUALITY(map.getNumGlobalEntries(), numGlobalEntries);
    TEST_EQUALITY_CONST(map.getNumMyEntries(), as<Ordinal>(2));
    TEST_EQUALITY_CONST(map.getIndexBase(), as<Ordinal>(indexBase));
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), as<Ordinal>(indexBase));
    TEST_EQUALITY_CONST(map.getMaxLocalIndex(), as<Ordinal>(1));
    TEST_EQUALITY_CONST(map.getMinGlobalIndex(), myGlobal[indexBase]);
    TEST_EQUALITY_CONST(map.getMaxGlobalIndex(), myGlobal[1]);
    TEST_EQUALITY_CONST(map.getMinAllGlobalIndex(), as<Ordinal>(indexBase));
    TEST_EQUALITY_CONST(map.getMaxAllGlobalIndex(), as<Ordinal>(numGlobalEntries-1));
    TEST_EQUALITY( map.getLocalIndex(myGlobal[0]), myLocal[0] );
    TEST_EQUALITY( map.getLocalIndex(myGlobal[1]), myLocal[1] );
    TEST_EQUALITY( map.getGlobalIndex(myLocal[0]), myGlobal[0] );
    TEST_EQUALITY( map.getGlobalIndex(myLocal[1]), myGlobal[1] );
    TEST_THROW( map.getGlobalIndex(2), std::invalid_argument);
    TEST_THROW( map.getLocalIndex(numGlobalEntries),  std::invalid_argument);  // all procs fail, all throw
    if (numImages > 1) {
      TEST_THROW( map.getLocalIndex(numGlobalEntries-1),std::invalid_argument);                   // all but one fail, all throw
    }
    TEST_COMPARE_ARRAYS( map.getMyGlobalEntries(), myGlobal);
    TEST_EQUALITY_CONST( map.isMyLocalIndex(0), true );
    TEST_EQUALITY_CONST( map.isMyLocalIndex(1), true );
    TEST_EQUALITY_CONST( map.isMyLocalIndex(2), false ); // just try a couple
    TEST_EQUALITY_CONST( map.isMyLocalIndex(3), false );
    for (Ordinal i=as<Ordinal>(0); i < numGlobalEntries; ++i) {
      if (find(myGlobal.begin(),myGlobal.end(),i) == myGlobal.end()) {
        TEST_EQUALITY_CONST( map.isMyGlobalIndex(i), false );
      }
      else {
        TEST_EQUALITY_CONST( map.isMyGlobalIndex(i), true );
      }
    }
    /* FINISH: Methods left to test
       void 	getRemoteIndexList (const std::vector< Ordinal > &GIDList, std::vector< Ordinal > &imageIDList, std::vector< Ordinal > &LIDList) const
       void 	getRemoteIndexList (const std::vector< Ordinal > &GIDList, std::vector< Ordinal > &imageIDList) const
     */
  }


  // 
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD


# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, invalidConstructor1, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, invalidConstructor2, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, invalidConstructor3, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, compatabilityTests, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, sameasTests, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, ContigUniformMap, ORDINAL )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, invalidConstructor1, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, invalidConstructor2, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, invalidConstructor3, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, compatabilityTests, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, sameasTests, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Map, ContigUniformMap, ORDINAL )

    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(LongInt)
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
