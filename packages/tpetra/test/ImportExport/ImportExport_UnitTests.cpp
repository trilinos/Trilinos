#include "Teuchos_UnitTestHarness.hpp"

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
// #include "Tpetra_Export.hpp" NOT FINISHED
#include "Tpetra_MultiVector.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::DefaultPlatform;
  using Tpetra::Platform;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Tpetra::Map;
  using Tpetra::Import;
  // using Tpetra::Export;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
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

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ImportExport, basic, Ordinal )
  {
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    // create a platform  
    RCP<const Platform<Ordinal> > platform = getDefaultPlatform<Ordinal>();
    // create a comm  
    //RCP<Comm<Ordinal> > comm = platform.createComm();
    //const int numImages = comm->getSize();
    //const int myImageID = comm->getRank();
    // create Maps
    Map<Ordinal> source(as<Ordinal>(-1),as<Ordinal>(10),ZERO,*platform),
                 target(as<Ordinal>(-1),as<Ordinal>(5) ,ZERO,*platform);
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
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ImportExport, GetNeighbors, Ordinal, Scalar )
  {
    typedef ScalarTraits<Scalar> ST;
    typedef Tpetra::MultiVector<Ordinal,Scalar> MV;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    const Ordinal ZERO = OrdinalTraits<Ordinal>::zero();
    const Ordinal ONE = OrdinalTraits<Ordinal>::one();
    const Ordinal NEGONE = ZERO - ONE;
    // create a platform  
    const Platform<Ordinal> & platform = *(getDefaultPlatform<Ordinal>());
    // create a comm  
    RCP<Comm<Ordinal> > comm = platform.createComm();
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
    Map<Ordinal> map(NEGONE,numLocal,indexBase,platform), 
                nmap(NEGONE,neighbors(),indexBase,platform);
    // mv = [ 0    n     2n    3n    4n ]
    //      [...  ....  ....  ....  ....]
    //      [n-1  2n-1  3n-1  4n-1  5n-1]
    MV mvMine(map,numVectors); 
    for (int j=0; j<numVectors; ++j) {
      mvMine.replaceMyValue(0,j,as<Scalar>(myImageID + j*numImages));
    }
    Import<Ordinal> importer(map,nmap);
    //Export<Ordinal> exporter(map,nmap);
    bool local_success = true;
    for (int t=0; t<1; t++) {
      // t=0 : use importer+doImport()
      // t=1 : use importer+doExport()
      // t=2 : use exporter+doImport()
      // t=3 : use exporter+doExport()
      MV mvWithNeighbors(nmap,numVectors);
      switch (t) {
        case 0:
          mvWithNeighbors.doImport(mvMine,importer,REPLACE);
          break;
        //case 1:
        //  mvMine.doExport(mvWithNeighbors,importer,REPLACE);
        //  break;
        //case 2:
        //  mvWithNeighbors.doImport(mvMine,exporter,REPLACE);
        //  break;
        //case 3:
        //  mvMine.doExport(mvWithNeighbors,exporter,REPLACE);
        //  break;
      }
      // test the numbers
      if (myImageID == 0) {
        for (int j=0; j<numVectors; ++j) {
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],0,as<Scalar>(j*numImages));           // me 
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
          TEST_ARRAY_ELE_EQUALITY(mvWithNeighbors[j],2,as<Scalar>(myImageID+j*numImages)-ST::one()); // neighbor
        }
      }
    }
    success &= local_success;
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
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ImportExport, GetNeighbors     , ORDINAL, SCALAR )

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

