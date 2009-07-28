#include "Teuchos_UnitTestHarness.hpp"

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::DefaultPlatform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Tpetra::Map;
  using Tpetra::Import;
  using Tpetra::Export;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
  using Tpetra::ADD;
  using std::ostream_iterator;
  using std::endl;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

#define PRINT_VECTOR(v) \
   { \
     out << #v << ": "; \
     copy(v.begin(), v.end(), ostream_iterator<Ordinal>(out," ")); \
     out << endl; \
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
      DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  // 

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ImportExport, basic, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    //const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();
    // create Maps
    Map<Ordinal> source(as<Ordinal>(-1),as<Ordinal>(10),ZERO,comm),
                 target(as<Ordinal>(-1),as<Ordinal>(5) ,ZERO,comm);
    // create Import object
    Import<Ordinal> importer(source, target);
    
    Ordinal same = importer.getNumSameIDs();
    Ordinal permute = importer.getNumPermuteIDs();
    Ordinal remote = importer.getNumRemoteIDs();
    Ordinal sum = same + permute + remote;
    Ordinal expectedSum = target.getNumMyEntries();
    TEST_EQUALITY( sum, expectedSum );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ImportExport, GetNeighborsForward, Ordinal, Scalar )
  {
    // import with the importer to duplicate
    // export with the exporter to add and reduce
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,Ordinal> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal  = ONE;
    const Ordinal numVectors = as<Ordinal>(5);
    // my neighbors: myImageID-1, me, myImageID+1
    Array<Ordinal> neighbors;
    if (myImageID != ZERO       ) neighbors.push_back(myImageID-1);
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) neighbors.push_back(myImageID+1);
    // two maps: one has one entries per node, the other is the 1-D neighbors
    Map<Ordinal> smap(NEGONE,numLocal,indexBase,comm), 
                 tmap(NEGONE,neighbors(),indexBase,comm);
    // mvMine = [myImageID  myImageID+numImages ... myImageID+4*numImages]
    MV mvMine(smap,numVectors); 
    for (int j=0; j<numVectors; ++j) {
      mvMine.replaceMyValue(0,j,as<Scalar>(myImageID + j*numImages));
    }
    // create Import from smap to tmap, Export from tmap to smap, test them
    Import<Ordinal> importer(smap,tmap);
    Export<Ordinal> exporter(tmap,smap);
    bool local_success = true;
    // importer testing: FINISH
    TEST_EQUALITY_CONST( importer.getSourceMap() == smap, true );
    TEST_EQUALITY_CONST( importer.getTargetMap() == tmap, true );
    TEST_EQUALITY( importer.getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
    TEST_EQUALITY( importer.getNumPermuteIDs(), (myImageID == 0 ? 0 : 1) );
    TEST_EQUALITY( importer.getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
    TEST_EQUALITY( importer.getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
    // exporter testing: FINISH
    TEST_EQUALITY_CONST( exporter.getSourceMap() == tmap, true );
    TEST_EQUALITY_CONST( exporter.getTargetMap() == smap, true );
    TEST_EQUALITY( importer.getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
    TEST_EQUALITY( exporter.getNumPermuteIDs(), (myImageID == 0 ? 0 : 1) );
    // import neighbors, test their proper arrival
    //                   [ 0    n     2n    3n    4n ]
    // mvWithNeighbors = [...  ....  ....  ....  ....]
    //                   [n-1  2n-1  3n-1  4n-1  5n-1]
    MV mvWithNeighbors(tmap,numVectors);
    mvWithNeighbors.doImport(mvMine,importer,REPLACE);
    if (myImageID == 0) {
      for (int j=0; j<numVectors; ++j) {
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(myImageID+j*numImages)); // me
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],1,as<Scalar>(j*numImages)+ST::one()); // neighbor
      }
    }
    else if (myImageID == numImages-1) {
      for (int j=0; j<numVectors; ++j) {
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],1,as<Scalar>(myImageID+j*numImages));           // me
      }
    }
    else {
      for (int j=0; j<numVectors; ++j) {
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],1,as<Scalar>(myImageID+j*numImages));           // me
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],2,as<Scalar>(myImageID+j*numImages)+ST::one()); // neighbor
      }
    }
    // export values, test 
    mvMine.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    mvMine.doExport(mvWithNeighbors,exporter,ADD);
    if (myImageID == 0 || myImageID == numImages-1) {
      for (int j=0; j<numVectors; ++j) {
        // contribution from me and one neighbor: double original value
        TEST_EQUALITY(mvMine[j][0],Teuchos::as<Scalar>(2.0)*as<Scalar>(myImageID+j*numImages));
      }
    }
    else {
      for (int j=0; j<numVectors; ++j) {
        // contribution from me and two neighbors: triple original value
        TEST_EQUALITY(mvMine[j][0],Teuchos::as<Scalar>(3.0)*as<Scalar>(myImageID+j*numImages));
      }
    }
    // 
    success &= local_success;
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ImportExport, GetNeighborsBackward, Ordinal, Scalar )
  {
    // import with the exporter to duplicate
    // export with the importer to add and reduce
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Scalar,Ordinal> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    // create a comm  
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    const Ordinal indexBase = ZERO;
    const Ordinal numLocal  = ONE;
    const Ordinal numVectors = as<Ordinal>(5);
    // my neighbors: myImageID-1, me, myImageID+1
    Array<Ordinal> neighbors;
    if (myImageID != ZERO       ) neighbors.push_back(myImageID-1);
    neighbors.push_back(myImageID);
    if (myImageID != numImages-1) neighbors.push_back(myImageID+1);
    // two maps: one has one entries per node, the other is the 1-D neighbors
    Map<Ordinal> smap(NEGONE,numLocal,indexBase,comm), 
                 tmap(NEGONE,neighbors(),indexBase,comm);
    // mvMine = [myImageID  myImageID+numImages ... myImageID+4*numImages]
    MV mvMine(smap,numVectors); 
    for (int j=0; j<numVectors; ++j) {
      mvMine.replaceMyValue(0,j,as<Scalar>(myImageID + j*numImages));
    }
    // create Import from smap to tmap, Export from tmap to smap, test them
    Import<Ordinal> importer(smap,tmap);
    Export<Ordinal> exporter(tmap,smap);
    bool local_success = true;
    // importer testing: FINISH
    TEST_EQUALITY_CONST( importer.getSourceMap() == smap, true );
    TEST_EQUALITY_CONST( importer.getTargetMap() == tmap, true );
    TEST_EQUALITY( importer.getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
    TEST_EQUALITY( importer.getNumPermuteIDs(), (myImageID == 0 ? 0 : 1) );
    TEST_EQUALITY( importer.getNumExportIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
    TEST_EQUALITY( importer.getNumRemoteIDs(), (myImageID == 0 || myImageID == numImages - 1 ? 1 : 2) );
    // exporter testing: FINISH
    TEST_EQUALITY_CONST( exporter.getSourceMap() == tmap, true );
    TEST_EQUALITY_CONST( exporter.getTargetMap() == smap, true );
    TEST_EQUALITY( importer.getNumSameIDs(), (myImageID == 0 ? 1 : 0) );
    TEST_EQUALITY( exporter.getNumPermuteIDs(), (myImageID == 0 ? 0 : 1) );
    // import neighbors, test their proper arrival
    //                   [ 0    n     2n    3n    4n ]
    // mvWithNeighbors = [...  ....  ....  ....  ....]
    //                   [n-1  2n-1  3n-1  4n-1  5n-1]
    MV mvWithNeighbors(tmap,numVectors);
    mvWithNeighbors.doImport(mvMine,exporter,REPLACE);
    if (myImageID == 0) {
      for (int j=0; j<numVectors; ++j) {
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(myImageID+j*numImages)); // me
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],1,as<Scalar>(j*numImages)+ST::one()); // neighbor
      }
    }
    else if (myImageID == numImages-1) {
      for (int j=0; j<numVectors; ++j) {
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],1,as<Scalar>(myImageID+j*numImages));           // me
      }
    }
    else {
      for (int j=0; j<numVectors; ++j) {
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],1,as<Scalar>(myImageID+j*numImages));           // me
        TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],2,as<Scalar>(myImageID+j*numImages)+ST::one()); // neighbor
      }
    }
    // export values, test 
    mvMine.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
    mvMine.doExport(mvWithNeighbors,importer,ADD);
    if (myImageID == 0 || myImageID == numImages-1) {
      for (int j=0; j<numVectors; ++j) {
        // contribution from me and one neighbor: double original value
        TEST_EQUALITY(mvMine[j][0],Teuchos::as<Scalar>(2.0)*as<Scalar>(myImageID+j*numImages));
      }
    }
    else {
      for (int j=0; j<numVectors; ++j) {
        // contribution from me and two neighbors: triple original value
        TEST_EQUALITY(mvMine[j][0],Teuchos::as<Scalar>(3.0)*as<Scalar>(myImageID+j*numImages));
      }
    }
    // 
    success &= local_success;
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, &globalSuccess_int );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  //
  // INSTANTIATIONS
  //

#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)\
     typedef std::complex<float> ComplexFloat; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexFloat)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)\
     typedef std::complex<double> ComplexDouble; \
     UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, ComplexDouble)
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ImportExport, GetNeighborsForward,  ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ImportExport, GetNeighborsBackward, ORDINAL, SCALAR )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ImportExport, basic, ORDINAL ) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, double) \
      UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL) 

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ImportExport, basic, ORDINAL ) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, float) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, double) \
      UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(ORDINAL) \
      UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(ORDINAL)

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

