#include "Teuchos_UnitTestHarness.hpp"

#include <map>
#include <vector>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
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
  using Tpetra::CrsMatrix;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
  using Tpetra::ADD;
  using std::ostream_iterator;
  using std::endl;
  using Tpetra::createContigMap;

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

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( CrsGraphImportExport, doImport, Ordinal ) 
  {
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm and node
    const RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize(),
              myImageID = comm->getRank();
    if (numImages < 2) return;

    //Create a Map that is evenly-distributed, and another that has all
    //elements on proc 0.
    {
      const int tgt_num_local_elements = 3;
      const int src_num_local_elements = (myImageID == 0 ? numImages*tgt_num_local_elements : 0);

      // create Maps
      RCP<const Map<Ordinal> >    src_map = createContigMap<Ordinal,Ordinal>(INVALID,src_num_local_elements,comm), 
                                  tgt_map = createContigMap<Ordinal,Ordinal>(INVALID,tgt_num_local_elements,comm);

      // create CrsGraph objects
      RCP<CrsGraph<Ordinal> >    src_graph = rcp(new CrsGraph<Ordinal>(src_map, 1)),
                                 tgt_graph = rcp(new CrsGraph<Ordinal>(tgt_map, 1));

      //Create a simple diagonal src-graph:
      Array<Ordinal> diag(1);
      Ordinal row = 0;
      for (size_t i=0; i<src_map->getNodeNumElements(); ++i, ++row) {
        Ordinal globalrow = src_map->getGlobalElement(row);
        diag[0] = globalrow;
        src_graph->insertGlobalIndices( globalrow, diag() );
      }

      Import<Ordinal> importer(src_map, tgt_map);
      tgt_graph->doImport(*src_graph, importer, Tpetra::INSERT);

      tgt_graph->fillComplete();

      //Now loop through tgt_graph and make sure it is diagonal:
      row = 0;
      for (size_t i=0; i<tgt_map->getNodeNumElements(); ++i, ++row) {
        ArrayRCP<const Ordinal> rowview = tgt_graph->getLocalRowView(row);
        TEST_EQUALITY(rowview.size(), 1);
        TEST_EQUALITY(rowview[0], row);
      }
    }


    //For the next test, we need an even number of processes:
    if (numImages%2 == 0) 
    {
      //Create Maps that are distributed differently but have the same global
      //number of elements. The source-map will have 3 elements on even-numbered
      //processors and 5 on odd-numbered processors. The target-map will have
      //4 elements on each processor:
      Ordinal src_num_local = 5;
      if (myImageID%2 == 0) src_num_local = 3;
      Ordinal tgt_num_local = 4;

      RCP<const Map<Ordinal> >    src_map = createContigMap<Ordinal,Ordinal>(INVALID,src_num_local,comm),  
                                  tgt_map = createContigMap<Ordinal,Ordinal>(INVALID,tgt_num_local,comm);  

      RCP<CrsGraph<Ordinal> >    src_graph = rcp(new CrsGraph<Ordinal>(src_map, 24)),
                                 tgt_graph = rcp(new CrsGraph<Ordinal>(tgt_map, 24));

      //This time make src_graph be a full lower-triangle:
      //Each row of column-indices will have length 'globalrow'+1, and contain
      //column-indices 0 .. 'globalrow'
      Array<Ordinal> cols(1);
      for (Ordinal globalrow=src_map->getMinGlobalIndex();
          globalrow<=src_map->getMaxGlobalIndex(); ++globalrow)
      {
        cols.resize(globalrow+1);
        for (Ordinal col=0; col<globalrow+1; ++col) {
          cols[col] = col;
        }
        src_graph->insertGlobalIndices( globalrow, cols() );
      }

      Import<Ordinal> importer(src_map, tgt_map);
      tgt_graph->doImport(*src_graph, importer, Tpetra::INSERT);

      src_graph->fillComplete();
      tgt_graph->fillComplete();

      //now we're going to loop through tgt_graph and make sure that
      //each row has length 'globalrow'+1 and has the correct contents:
      const Teuchos::RCP<const Map<Ordinal> >& colmap = tgt_graph->getColMap();

      for (Ordinal globalrow=tgt_map->getMinGlobalIndex(); globalrow<=tgt_map->getMaxGlobalIndex(); ++globalrow)
      {
        Ordinal localrow = tgt_map->getLocalElement(globalrow);
        ArrayRCP<const Ordinal> rowview = tgt_graph->getLocalRowView(localrow);
        TEST_EQUALITY(rowview.size(), globalrow+1);
        for (Ordinal j=0; j<globalrow+1; ++j) {
          TEST_EQUALITY(colmap->getGlobalElement(rowview[j]), j);
        }
      }
    }
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrixImportExport, doImport, Ordinal, Scalar ) 
  {
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize(),
              myImageID = comm->getRank();
    // need at least two processes to test import/export
    if (numImages < 2) return;

    {
      // Create a Map that is evenly-distributed, and another that has all
      // elements on proc 0.
      const Ordinal tgt_num_local_elements = 3;
      const Ordinal src_num_local_elements = ( myImageID == 0 ? numImages*tgt_num_local_elements : 0 );

      // create Maps
      RCP<const Map<Ordinal> > src_map = createContigMap<Ordinal,Ordinal>(INVALID,src_num_local_elements,comm), 
                               tgt_map = createContigMap<Ordinal,Ordinal>(INVALID,tgt_num_local_elements,comm);

      // create CrsGraph objects
      RCP<CrsMatrix<Scalar,Ordinal> > src_mat = rcp(new CrsMatrix<Scalar,Ordinal>(src_map, 1)),
                                      tgt_mat = rcp(new CrsMatrix<Scalar,Ordinal>(tgt_map, 1));

      // Create a simple diagonal src-graph:
      for (Ordinal globalrow = src_map->getMinGlobalIndex(); globalrow <= src_map->getMaxGlobalIndex(); ++globalrow) 
      {
        src_mat->insertGlobalValues( globalrow, tuple<Ordinal>(globalrow), tuple<Scalar>(globalrow) );
      }

      // Create the importer
      Import<Ordinal> importer(src_map, tgt_map);
      // Do the import and complete fill on tgt_map
      tgt_mat->doImport(*src_mat, importer, Tpetra::INSERT);
      tgt_mat->fillComplete();

      // Now loop through tgt_mat and make sure it is diagonal:
      for (Ordinal localrow = tgt_map->getMinLocalIndex(); localrow <= tgt_map->getMaxLocalIndex(); ++localrow) 
      {
        ArrayRCP<const Ordinal> rowinds;
        ArrayRCP<const Scalar>  rowvals;
        tgt_mat->getLocalRowView(localrow, rowinds, rowvals);
        TEST_EQUALITY_CONST(rowinds.size(), 1);
        TEST_EQUALITY(rowinds[0], as<Ordinal>(localrow));
        TEST_EQUALITY_CONST(rowvals.size(), 1);
        TEST_EQUALITY(rowvals[0], as<Scalar>(localrow));
      }
    }


    // For the next test, we need an even number of processes:
    if (numImages%2 == 0) {
      // Create Maps that are distributed differently but have the same global
      // number of elements. The source-map will have 3 elements on even-numbered
      // processors and 5 on odd-numbered processors. The target-map will have
      // 4 elements on each processor:
      const Ordinal src_num_local = (myImageID%2 == 0 ? 3 : 5);
      const Ordinal tgt_num_local = 4;

      RCP<const Map<Ordinal> > src_map = createContigMap<Ordinal,Ordinal>(INVALID,src_num_local,comm),  
                               tgt_map = createContigMap<Ordinal,Ordinal>(INVALID,tgt_num_local,comm);  

      RCP<CrsMatrix<Scalar,Ordinal> > src_mat = rcp(new CrsMatrix<Scalar,Ordinal>(src_map, 24)),
                                      tgt_mat = rcp(new CrsMatrix<Scalar,Ordinal>(tgt_map, 24));

      // This time make src_mat be a full lower-triangle:
      // Each row of column-indices will have length 'globalrow', and contain
      // column-indices 0 .. 'globalrow'-1
      Array<Ordinal> cols(1);
      Array<Scalar>  vals(1);
      for (Ordinal globalrow=src_map->getMinGlobalIndex(); globalrow<=src_map->getMaxGlobalIndex(); ++globalrow)
      {
        if (globalrow > 0) {
          cols.resize(globalrow);
          vals.resize(globalrow);
          for (Ordinal col=0; col<globalrow; ++col) {
            cols[col] = as<Ordinal>(col);
            vals[col] = as<Scalar>(col);
          }
          src_mat->insertGlobalValues( globalrow, cols(), vals() );
        }
      }

      Import<Ordinal> importer(src_map, tgt_map);
      tgt_mat->doImport(*src_mat, importer, Tpetra::INSERT);

      src_mat->fillComplete();
      tgt_mat->fillComplete();

      // now we're going to loop through tgt_mat and make sure that
      // each row has length 'globalrow' and has the correct contents:
      const Teuchos::RCP<const Map<Ordinal> > colmap = tgt_mat->getColMap();

      for (Ordinal globalrow=tgt_map->getMinGlobalIndex(); globalrow<=tgt_map->getMaxGlobalIndex(); ++globalrow)
      {
        Ordinal localrow = tgt_map->getLocalElement(globalrow);
        ArrayRCP<const Ordinal> rowinds;
        ArrayRCP<const Scalar> rowvals;
        tgt_mat->getLocalRowView(localrow, rowinds, rowvals);
        TEST_EQUALITY(rowinds.size(), globalrow);
        TEST_EQUALITY(rowvals.size(), globalrow);
        for (Teuchos_Ordinal j=0; j<rowinds.size(); ++j) {
          TEST_EQUALITY( colmap->getGlobalElement(rowinds[j]), as<Ordinal>(j) );
          TEST_EQUALITY( rowvals[j], as<Scalar>(j)  );
        }
      }
    }
  }



  //
  // INSTANTIATIONS
  //

#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( CrsGraphImportExport, doImport, ORDINAL )

#   define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrixImportExport, doImport, ORDINAL, SCALAR )

  // test CrsGraph for <int,int,DefaultNode>
  // if explicit instantiation is enabled, this configuration is always built
  // if not, it is implicitly instantiated
  // therefore, this is always possible
  UNIT_TEST_GROUP_ORDINAL(int)

  // test CrsMatrix for some scalar
  // if explicit instantiation is enabled, test for all enabled types
  // if not, test for float
#if defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION)
  #if defined(HAVE_TPETRA_INST_DOUBLE)
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int,double)
  #elif defined(HAVE_TPETRA_INST_FLOAT)
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int,float)
  #elif defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
    typedef std::complex<float> ComplexFloat; \
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int, ComplexFloat)
  #elif defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
    typedef std::complex<double> ComplexDouble; \
    UNIT_TEST_GROUP_ORDINAL_SCALAR(int, ComplexDouble)
  #endif
#else
  UNIT_TEST_GROUP_ORDINAL_SCALAR(int,float)
#endif

}

