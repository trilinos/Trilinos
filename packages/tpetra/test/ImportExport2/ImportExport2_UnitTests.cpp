#include "Teuchos_UnitTestHarness.hpp"

#include <map>
#include <vector>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_CrsGraph.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Tpetra::DefaultPlatform;
  using Tpetra::global_size_t;
  using Teuchos::arrayViewFromVector;
  using Teuchos::OrdinalTraits;
  using Teuchos::tuple;
  using Tpetra::Map;
  using Tpetra::Import;
  using Tpetra::Export;
  using Tpetra::CrsGraph;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
  using Tpetra::ADD;
  using std::ostream_iterator;
  using std::endl;

  typedef DefaultPlatform::DefaultPlatformType::NodeType Node;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

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

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraphImportExport, doImport, Ordinal ) {
    const Ordinal indexBase = OrdinalTraits<Ordinal>::zero();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize(),
              myImageID = comm->getRank();
    if (numImages < 2) return;

    //Create a Map that is evenly-distributed, and another that has all
    //elements on proc 0.

    Ordinal target_num_local_elements = 3;
    Ordinal src_num_local_elements = 0;
    if (myImageID == 0) src_num_local_elements = numImages*target_num_local_elements;

    // create Maps
    RCP<const Map<Ordinal,Ordinal,Node> > src_map = rcp(new Map<Ordinal,Ordinal,Node>(INVALID,src_num_local_elements,indexBase,comm) ),
                                          target_map = rcp(new Map<Ordinal,Ordinal,Node>(INVALID, target_num_local_elements,indexBase,comm) );

    // create CrsGraph objects
    RCP<CrsGraph<Ordinal,Ordinal,Node> > src_graph = rcp(new CrsGraph<Ordinal,Ordinal,Node>(src_map, 1));
    RCP<CrsGraph<Ordinal,Ordinal,Node> > target_graph = rcp(new CrsGraph<Ordinal,Ordinal,Node>(target_map, 1));

    //Create a simple diagonal src-graph:
    Teuchos::Array<Ordinal> diag(1);
    Ordinal row = 0;
    for(size_t i=0; i<src_map->getNodeNumElements(); ++i, ++row) {
      Ordinal globalrow = src_map->getGlobalElement(row);
      diag[0] = globalrow;
      src_graph->insertGlobalIndices( globalrow, diag() );
    }

    Import<Ordinal> importer(src_map, target_map);
    target_graph->doImport(*src_graph, importer, Tpetra::INSERT);

    target_graph->fillComplete();

    //Now loop through target_graph and make sure it is diagonal:
    row = 0;
    for(size_t i=0; i<target_map->getNodeNumElements(); ++i, ++row) {
      Teuchos::ArrayRCP<const Ordinal> rowview = target_graph->getLocalRowView(row);
      TEST_EQUALITY(rowview.size(), 1);
      TEST_EQUALITY(rowview[0], row);
    }


    //For the next test, we need an even number of processes:
    if (numImages%2 != 0) return;

    //Create Maps that are distributed differently but have the same global
    //number of elements. The source-map will have 3 elements on even-numbered
    //processors and 5 on odd-numbered processors. The target-map will have
    //4 elements on each processor:
    Ordinal src_num_local = 5;
    if (myImageID%2 == 0) src_num_local = 3;
    Ordinal target_num_local = 4;

    RCP<const Map<Ordinal,Ordinal,Node> > src_map2 = rcp(new Map<Ordinal,Ordinal,Node>(INVALID,src_num_local,indexBase,comm) ),
                                          target_map2 = rcp(new Map<Ordinal,Ordinal,Node>(INVALID,target_num_local,indexBase,comm) );

    RCP<CrsGraph<Ordinal,Ordinal,Node> > src_graph2 = rcp(new CrsGraph<Ordinal,Ordinal,Node>(src_map2, 24));
    RCP<CrsGraph<Ordinal,Ordinal,Node> > target_graph2 = rcp(new CrsGraph<Ordinal,Ordinal,Node>(target_map2, 24));

    //This time make src_graph2 be a full lower-triangle:
    //Each row of column-indices will have length 'globalrow'+1, and contain
    //column-indices 0 .. 'globalrow'
    Teuchos::Array<Ordinal> cols(1);
    for(Ordinal globalrow=src_map2->getMinGlobalIndex();
                globalrow<=src_map2->getMaxGlobalIndex(); ++globalrow)
    {
      cols.resize(globalrow+1);
      for(Ordinal col=0; col<globalrow+1; ++col) {
        cols[col] = col;
      }
      src_graph2->insertGlobalIndices( globalrow, cols() );
    }

    Import<Ordinal> importer2(src_map2, target_map2);
    target_graph2->doImport(*src_graph2, importer2, Tpetra::INSERT);

    src_graph2->fillComplete();
    target_graph2->fillComplete();

    //now we're going to loop through target_graph2 and make sure that
    //each row has length 'globalrow'+1 and has the correct contents:
    const Teuchos::RCP<const Map<Ordinal,Ordinal,Node> >& colmap =
       target_graph2->getColMap();

    for(Ordinal globalrow=target_map2->getMinGlobalIndex();
                globalrow<=target_map2->getMaxGlobalIndex(); ++globalrow)
    {
      Ordinal localrow = target_map2->getLocalElement(globalrow);
      Teuchos::ArrayRCP<const Ordinal> rowview = target_graph2->getLocalRowView(localrow);
      TEST_EQUALITY(rowview.size(), globalrow+1);
      for(Ordinal j=0; j<globalrow+1; ++j) {
        TEST_EQUALITY(colmap->getGlobalElement(rowview[j]), j);
      }
    }
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

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraphImportExport, doImport, ORDINAL )

    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraphImportExport, doImport, ORDINAL )

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

