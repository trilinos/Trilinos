#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixSolveOp.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "TpetraExt_MatrixMatrix.hpp"

#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOS_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  using std::endl;
  using std::swap;

  using std::string;

  using Teuchos::TypeTraits::is_same;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  using Tpetra::createCrsMatrixSolveOp;
  using Tpetra::createCrsMatrixMultiplyOp;
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  using Kokkos::SerialNode;
  RCP<SerialNode> snode;
#ifdef HAVE_KOKKOS_TBB
  using Kokkos::TBBNode;
  RCP<TBBNode> tbbnode;
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
  using Kokkos::TPINode;
  RCP<TPINode> tpinode;
#endif
#ifdef HAVE_KOKKOS_THRUST
  using Kokkos::ThrustGPUNode;
  RCP<ThrustGPUNode> thrustnode;
#endif

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  string filedir;

#define DYN_ARRAYVIEW_TO_ARRAY(Type, arr, av) \
  { \
    ArrayView<Type> av2 = av; \
    arr.resize(av2.size()); \
    arr.assign(av2.begin(),av2.end()); \
  }

#define STD_TESTS(matrix) \
  { \
    using Teuchos::outArg; \
    RCP<const Comm<int> > STCOMM = matrix.getComm(); \
    ArrayView<const GO> STMYGIDS = matrix.getRowMap()->getNodeElementList(); \
    ArrayView<const LO> loview; \
    ArrayView<const Scalar> sview; \
    size_t STMAX = 0; \
    for (size_t STR=0; STR < matrix.getNodeNumRows(); ++STR) { \
      const size_t numEntries = matrix.getNumEntriesInLocalRow(STR); \
      TEST_EQUALITY( numEntries, matrix.getNumEntriesInGlobalRow( STMYGIDS[STR] ) ); \
      matrix.getLocalRowView(STR,loview,sview); \
      TEST_EQUALITY( static_cast<size_t>(loview.size()), numEntries ); \
      TEST_EQUALITY( static_cast<size_t>( sview.size()), numEntries ); \
      STMAX = std::max( STMAX, numEntries ); \
    } \
    TEST_EQUALITY( matrix.getNodeMaxNumRowEntries(), STMAX ); \
    global_size_t STGMAX; \
    Teuchos::reduceAll<int,global_size_t>( *STCOMM, Teuchos::REDUCE_MAX, STMAX, outArg(STGMAX) ); \
    TEST_EQUALITY( matrix.getGlobalMaxNumRowEntries(), STGMAX ); \
  }


  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption(
        "filedir",&filedir,"Directory of expected matrix files.");
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
    RCP<const Comm<int> > ret;
    if (testMpi) {
      ret = DefaultPlatform::getDefaultPlatform().getComm();
    }
    else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  template <class Node>
  RCP<Node> getNode() {
    assert(false);
  }

  template <>
  RCP<SerialNode> getNode<SerialNode>() {
    if (snode == null) {
      Teuchos::ParameterList pl;
      snode = rcp(new SerialNode(pl));
    }
    return snode;
  }

#ifdef HAVE_KOKKOS_TBB
  template <>
  RCP<TBBNode> getNode<TBBNode>() {
    if (tbbnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tbbnode = rcp(new TBBNode(pl));
    }
    return tbbnode;
  }
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
  template <>
  RCP<TPINode> getNode<TPINode>() {
    if (tpinode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      tpinode = rcp(new TPINode(pl));
    }
    return tpinode;
  }
#endif

#ifdef HAVE_KOKKOS_THRUST
  template <>
  RCP<ThrustGPUNode> getNode<ThrustGPUNode>() {
    if (thrustnode == null) {
      Teuchos::ParameterList pl;
      pl.set<int>("Num Threads",0);
      pl.set<int>("Verbose",1);
      thrustnode = rcp(new ThrustGPUNode(pl));
    }
    return thrustnode;
  }
#endif

  //
  // UNIT TESTS
  // 

  #include "CrsMatrix_NonlocalAfterResume.hpp"

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, BadCalls, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef RCP<const Map<LO,GO,Node> > RCPMap;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    // create the zero matrix
    RCP<CrsMatrix<Scalar,LO,GO,Node> > zero;
    {
      RCPMap map  = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
      MV mv(map,1);
      zero = rcp( new MAT(map,0,DynamicProfile) );
      TEST_THROW(zero->apply(mv,mv), std::runtime_error);
#   if defined(HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS)
      // throw exception because we required increased allocation
      TEST_THROW(zero->insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(0),tuple<Scalar>(ST::one())), std::runtime_error);
#   endif
      TEST_EQUALITY_CONST( zero->getProfileType() == DynamicProfile, true );
      zero->fillComplete();
    }
    STD_TESTS((*zero));
    TEST_EQUALITY_CONST( zero->getRangeMap() == zero->getDomainMap(), true );
    TEST_EQUALITY_CONST( zero->getFrobeniusNorm(), MT::zero() );
    const RCPMap drmap = zero->getDomainMap();
    {
      MV mv1(drmap,1), mv2(drmap,2), mv3(drmap,3);
      TEST_THROW(zero->apply(mv2,mv1), std::runtime_error); // MVs have different number of vectors
      TEST_THROW(zero->apply(mv2,mv3), std::runtime_error); // MVs have different number of vectors
    }
    // test that our assumptions on the maps are correct: 
    // that is, that badmap is not equal to the range, domain, row or colum map of the matrix
    const RCPMap badmap = createContigMapWithNode<LO,GO>(INVALID,1,comm,node);
    TEST_EQUALITY_CONST( badmap != zero->getRowMap(), true );
    TEST_EQUALITY_CONST( badmap != zero->getColMap(), true );
    TEST_EQUALITY_CONST( badmap != zero->getDomainMap(), true );
    TEST_EQUALITY_CONST( badmap != zero->getRangeMap(),  true );
    TEST_EQUALITY_CONST( *badmap != *zero->getRowMap(), true );
    TEST_EQUALITY_CONST( *badmap != *zero->getColMap(), true );
    TEST_EQUALITY_CONST( *badmap != *zero->getDomainMap(), true );
    TEST_EQUALITY_CONST( *badmap != *zero->getRangeMap(),  true );
    // now test the multivector against the matrix operators
    // Bugzilla bug #5247
    {
      MV mvbad(badmap,1);
#ifdef HAVE_TPETRA_DEBUG
      const Scalar ONE = ST::one(), ZERO = ST::zero();
      // tests in localSolve() and localMultiply() are only done in a debug build
      MV mvcol(zero->getColMap(),1);
      MV mvrow(zero->getRowMap(),1);
      TEST_THROW(zero->template localMultiply<Scalar>(mvcol,mvbad,  NO_TRANS,ONE,ZERO), std::runtime_error); // bad output map
      TEST_THROW(zero->template localMultiply<Scalar>(mvbad,mvrow,  NO_TRANS,ONE,ZERO), std::runtime_error); // bad input map
      TEST_THROW(zero->template localMultiply<Scalar>(mvbad,mvcol,CONJ_TRANS,ONE,ZERO), std::runtime_error); // bad output map
      TEST_THROW(zero->template localMultiply<Scalar>(mvrow,mvbad,CONJ_TRANS,ONE,ZERO), std::runtime_error); // bad input map
#endif
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EmptyFillComplete, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Vector<Scalar,LO,GO,Node> V;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      // create static-profile matrix, fill-complete without inserting (and therefore, without allocating)
      MAT matrix(map,1,StaticProfile);
      matrix.fillComplete(DoOptimizeStorage);
    }
    {
      // create dynamic-profile matrix, fill-complete without inserting (and therefore, without allocating)
      MAT matrix(map,1,DynamicProfile);
      matrix.fillComplete(DoOptimizeStorage);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, AdvancedGraphUsage, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Vector<Scalar,LO,GO,Node> V;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      CrsGraph<LO,GO,Node> diaggraph(map,1,StaticProfile);
      // A pre-constructed graph must be fill complete before being used to construct a CrsMatrix
      TEST_THROW( MAT matrix(rcpFromRef(diaggraph)), std::runtime_error );
    }
    {
      // create a simple diagonal graph
      CrsGraph<LO,GO,Node> diaggraph(map,1,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }
      // fill-complete the graph, but do not optimize the storage
      diaggraph.fillComplete(DoNotOptimizeStorage);
      TEST_EQUALITY_CONST( diaggraph.isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph.isStorageOptimized(), false );
      // matrix constructed with non-storage-optimized graph
      MAT mat1(rcpFromRef(diaggraph));
      // fill complete the matrix and ask it to optimize storage.
      // this is not allowed on a static graph, and will either throw an exception or ignore the request to optimize storage.
#ifdef HAVE_TPETRA_THROW_ABUSE_WARNINGS
      TEST_THROW( mat1.fillComplete(DoOptimizeStorage), std::runtime_error );
      TEST_EQUALITY_CONST( mat1.isFillComplete(), false );
      TEST_EQUALITY_CONST( mat1.isStorageOptimized(), false );
#else
      mat1.fillComplete(DoOptimizeStorage);
      TEST_EQUALITY_CONST( mat1.isFillComplete(), true );
      TEST_EQUALITY_CONST( mat1.isStorageOptimized(), false );
#endif
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithGraph, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef CrsGraph<LO,GO,Node>        GRPH;
    typedef Vector<Scalar,LO,GO,Node> V;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar SONE  = ST::one();
    const Scalar SZERO = ST::zero();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map with numLocal entries per node
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      //////////////////////////////////
      // create a simple tridiagonal graph
      GRPH trigraph(map,3,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        if (r == map->getMinAllGlobalIndex()) {
          trigraph.insertGlobalIndices(r,tuple(r,r+1));
        }
        else if (r == map->getMaxAllGlobalIndex()) {
          trigraph.insertGlobalIndices(r,tuple(r-1,r));
        }
        else {
          trigraph.insertGlobalIndices(r,tuple(r-1,r,r+1));
        }
      }
      trigraph.fillComplete(DoOptimizeStorage);
      // create a matrix using the tri-diagonal graph and test allowed functionality
      MAT matrix(rcpFromRef(trigraph));
      TEST_EQUALITY_CONST( matrix.getProfileType() == StaticProfile, true );
      // insert throws exception: not allowed with static graph
      TEST_THROW( matrix.insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(map->getMinGlobalIndex()),tuple(ST::one())), std::runtime_error );
      // suminto and replace are allowed
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        if (r == map->getMinAllGlobalIndex()) {
          matrix.replaceGlobalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
        }
        else if (r == map->getMaxAllGlobalIndex()) {
          matrix.replaceGlobalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
        }
        else {
          matrix.replaceGlobalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
        }
      }
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        // increment the diagonals
        matrix.sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
      }
      matrix.fillComplete();
      TEST_EQUALITY( matrix.getNodeNumDiags(), numLocal );
      TEST_EQUALITY( matrix.getGlobalNumDiags(), numImages*numLocal );
      TEST_EQUALITY( matrix.getGlobalNumEntries(), 3*numImages*numLocal - 2 );
      V dvec(map,false);
      dvec.randomize();
      matrix.getLocalDiagCopy(dvec);
      Array<Scalar> expectedDiags(numLocal, static_cast<Scalar>(2));
      ArrayRCP<const Scalar> dvec_view = dvec.get1dView();
      TEST_COMPARE_FLOATING_ARRAYS( expectedDiags(), dvec_view, MT::zero() );
    }
    {
      // create a simple diagonal graph
      GRPH diaggraph(map,1,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }
      diaggraph.fillComplete(DoOptimizeStorage);
      // Bug verification:
      //  Tpetra::CrsMatrix constructed with a graph was experiencing a seg-fault if setAllToScalar is called before 
      //  some other call allocates memory. This was because setAllToScalar has an incorrect if-statement that is 
      //  not allocating memory.
      // This bug has been fixed. Furthermore, CrsMatrix no longer utilizes lazy allocation when constructed with a graph. 
      // However, we will leave this test in place, because it still demonstrates valid behavior.
      MAT matrix(rcpFromRef(diaggraph));
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );
    }
    {
      // create a simple diagonal graph
      GRPH diaggraph(map,1,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph.insertGlobalIndices(r,tuple(r));
      }
      diaggraph.fillComplete(DoOptimizeStorage);
      TEST_EQUALITY_CONST( diaggraph.isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph.isStorageOptimized(), true );
      TEST_EQUALITY_CONST( diaggraph.isUpperTriangular(), true );
      TEST_EQUALITY_CONST( diaggraph.isLowerTriangular(), true );
      // Bug verification:
      // Tpetra::CrsMatrix constructed with a Optimized, Fill-Complete graph will not call fillLocalMatrix() 
      // in optimizeStorage(), because it returns early due to picking up the storage optimized bool from the graph.
      // As a result, the local mat-vec and mat-solve operations are never initialized, and localMultiply() and localSolve() 
      // fail with a complaint regarding the initialization of these objects.
      MAT matrix(rcpFromRef(diaggraph));
      TEST_NOTHROW( matrix.setAllToScalar( ST::one() ) );
      matrix.fillComplete(DoOptimizeStorage);
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
      TEST_EQUALITY_CONST( matrix.isUpperTriangular(), true );
      TEST_EQUALITY_CONST( matrix.isLowerTriangular(), true );
      // init x to ones(); multiply into y, solve in-situ in y, check result
      V x(map,false), y(map,false);
      x.putScalar(SONE);
      TEST_NOTHROW( matrix.localMultiply(x,y,NO_TRANS,SONE,SZERO) );
      TEST_NOTHROW( matrix.localSolve(y,y,NO_TRANS) );
      ArrayRCP<const Scalar> x_view = x.get1dView();
      ArrayRCP<const Scalar> y_view = y.get1dView();
      TEST_COMPARE_FLOATING_ARRAYS( y_view, x_view, MT::zero() );
    }
    {
      // create a simple diagonal graph
      RCP<GRPH> diaggraph = rcp( new GRPH(map,1,StaticProfile) );
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        diaggraph->insertGlobalIndices(r,tuple(r));
      }
      diaggraph->fillComplete(DoOptimizeStorage);
      TEST_EQUALITY_CONST( diaggraph->isFillComplete(), true );
      TEST_EQUALITY_CONST( diaggraph->isStorageOptimized(), true );
      TEST_EQUALITY_CONST( diaggraph->isUpperTriangular(), true );
      TEST_EQUALITY_CONST( diaggraph->isLowerTriangular(), true );
      // construct a matrix with the graph from another matrix
      MAT matrix1(diaggraph);
      TEST_EQUALITY( matrix1.getCrsGraph(), diaggraph );
      MAT matrix2( matrix1.getCrsGraph() );
      TEST_EQUALITY( matrix2.getCrsGraph(), matrix1.getCrsGraph() );
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithColMap, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar SONE  = ST::one();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map with numLocal entries per node using a pre-existing column map.
    // ensure:
    // * that the matrix uses this col map
    // * that it performs filtering during insertions
    // * that we can perform local or global insertions
    const size_t numLocal = 10; TEST_FOR_EXCEPTION( numLocal < 2, std::logic_error, "Test assumes that numLocal be greater than 1.");
    // these maps are equalivalent, but we should keep two distinct maps just to verify the general use case.
    RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    RCP<const Map<LO,GO,Node> > cmap = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    //////////////////////////////////
    // add tridiagonal entries, but use a diagonal column map.
    // result should be block diagonal matrix, with no importer/exporter.
    //
    // run this test twice; once where we insert global indices and once where we insert local indices
    // both are allowed with a specified column map; however, we can only test one at a time.
    // 
    // the first time, use const NNZ
    // the second, use NNZ array
    {
      MAT bdmat(rmap,cmap,3,StaticProfile);
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        // use global for the first one to verify that the matrix allows it
        // r-1 might be invalid, but the column map filtering should address that.
        bdmat.insertGlobalValues(r,tuple<GO>(r-1,r,r+1),tuple<Scalar>(SONE,SONE,SONE));
      }
      TEST_NOTHROW(bdmat.fillComplete());
      // nothing should have changed with regard to the row and column maps of the matrix
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      // check that filtering happened
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        if (r == rmap->getMinGlobalIndex() || r == rmap->getMaxGlobalIndex()) {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 2);
        }
        else {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 3);
        }
      }
    }
    {
      ArrayRCP<size_t> nnzperrow = arcp<size_t>(numLocal);
      std::fill(nnzperrow.begin(), nnzperrow.end(), 3);
      MAT bdmat(rmap,cmap,nnzperrow,StaticProfile);
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        // use local for the rest. need the column map
        // column map and row map are the same, so we only have to do one translation
        LO lid = cmap->getLocalElement(r);
        // as above, filtering via column map (required to happen for local and global) will save us for the invalid r-1
        bdmat.insertLocalValues(lid,tuple<LO>(lid-1,lid,lid+1),tuple<Scalar>(SONE,SONE,SONE));
      }
      TEST_NOTHROW(bdmat.fillComplete());
      // nothing should have changed with regard to the row and column maps of the matrix
      TEST_EQUALITY(bdmat.getRowMap(), rmap);
      TEST_EQUALITY_CONST(bdmat.hasColMap(), true);
      TEST_EQUALITY(bdmat.getColMap(), cmap);
      // check that filtering happened
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        if (r == rmap->getMinGlobalIndex() || r == rmap->getMaxGlobalIndex()) {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 2);
        }
        else {
          TEST_EQUALITY_CONST(bdmat.getNumEntriesInGlobalRow(r), 3);
        }
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, WithGraph_replaceLocal, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // generate a tridiagonal matrix
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Vector<Scalar,LO,GO,Node> V;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    CrsGraph<LO,GO,Node> graph(map,3,StaticProfile);
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      if (r == map->getMinAllGlobalIndex()) {
        graph.insertGlobalIndices(r,tuple(r,r+1));
      }
      else if (r == map->getMaxAllGlobalIndex()) {
        graph.insertGlobalIndices(r,tuple(r-1,r));
      }
      else {
        graph.insertGlobalIndices(r,tuple(r-1,r,r+1));
      }
    }
    graph.fillComplete(DoOptimizeStorage);
    // create a matrix using the graph
    MAT matrix(rcpFromRef(graph));
    TEST_EQUALITY_CONST( matrix.getProfileType() == StaticProfile, true );
    // insert throws exception: not allowed with static graph
    TEST_THROW( matrix.insertGlobalValues(map->getMinGlobalIndex(),tuple<GO>(map->getMinGlobalIndex()),tuple(ST::one())), std::runtime_error );
    // suminto and replace are allowed
    for (LO r=map->getMinLocalIndex(); r <= map->getMaxLocalIndex(); ++r) {
      if (r == map->getMinLocalIndex()) {
        matrix.replaceLocalValues(r, tuple(r,r+1), tuple(ST::one(),ST::one()) );
      }
      else if (r == map->getMaxLocalIndex()) {
        matrix.replaceLocalValues(r, tuple(r-1,r), tuple(ST::one(),ST::one()) );
      }
      else {
        matrix.replaceLocalValues(r, tuple(r-1,r,r+1), tuple(ST::one(),ST::one(),ST::one()) );
      }
    }
    for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
      // increment the diagonals
      matrix.sumIntoGlobalValues(r, tuple(r), tuple(ST::one()) );
    }
    matrix.fillComplete();
    TEST_EQUALITY( matrix.getNodeNumDiags(), numLocal );
    TEST_EQUALITY( matrix.getGlobalNumDiags(), numImages*numLocal );
    TEST_EQUALITY( matrix.getGlobalNumEntries(), 3*numImages*numLocal - 2 );
    V dvec(map,false);
    dvec.randomize();
    matrix.getLocalDiagCopy(dvec);
    Array<Scalar> expectedDiags(numLocal, static_cast<Scalar>(2));
    ArrayRCP<const Scalar> dvec_view = dvec.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS( expectedDiags(), dvec_view, MT::zero() );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ExceedStaticAlloc, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      MAT matrix(map,1,StaticProfile);
      // room for one on each row
      for (GO r=map->getMinGlobalIndex(); r<=map->getMaxGlobalIndex(); ++r) 
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // no room for any more
      GO r = map->getMinGlobalIndex();
      TEST_THROW( matrix.insertGlobalValues( r, tuple(r+1), tuple(ST::one()) ), std::runtime_error );
    }
    if (numImages > 1) {
      // add too many entries globally
      MAT matrix(map,1,StaticProfile);
      // room for one on each row
      for (GO r=map->getMinGlobalIndex(); r<=map->getMaxGlobalIndex(); ++r) 
      {
        matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
      }
      // always room for non-locals
      GO r = map->getMaxGlobalIndex() + 1;
      if (r > map->getMaxAllGlobalIndex()) r = map->getMinAllGlobalIndex();
      TEST_NOTHROW( matrix.insertGlobalValues( r, tuple(r), tuple(ST::one()) ) );
      // after communicating non-locals, failure trying to add them
      TEST_THROW( matrix.globalAssemble(), std::runtime_error );
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MultipleFillCompletes, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // test that an exception is thrown when we exceed statically allocated memory
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> ST;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    // create a Map
    const size_t numLocal = 1; // change to 10
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      // room for two on each row
      MAT matrix(map,2,StaticProfile);
      for (GO r=map->getMinGlobalIndex(); r <= map->getMaxGlobalIndex(); ++r) {
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
        TEST_NOTHROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())) );
      }
      // no room for more
      {
        GO r = map->getMinGlobalIndex();
        TEST_THROW( matrix.insertGlobalValues(r,tuple(r),tuple(ST::one())), std::runtime_error );
      }
      // fill complete adds them
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      TEST_NOTHROW( matrix.fillComplete(DoNotOptimizeStorage) );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
      // now there is room for more
      matrix.resumeFill();
      for (LO r=0; r<static_cast<LO>(numLocal); ++r) 
      {
        matrix.insertLocalValues(r,tuple(r),tuple(ST::one()));
      }
      TEST_NOTHROW( matrix.fillComplete(DoOptimizeStorage) );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
      // test that matrix is 3*I
      TEST_EQUALITY( matrix.getGlobalNumDiags(), numLocal*numImages );
      TEST_EQUALITY( matrix.getNodeNumDiags(), numLocal );
      TEST_EQUALITY( matrix.getGlobalNumEntries(), numLocal*numImages );
      TEST_EQUALITY( matrix.getNodeNumEntries(), numLocal );
      for (LO r=0; r<static_cast<LO>(numLocal); ++r) {
        ArrayView<const LO> inds;
        ArrayView<const Scalar> vals;
        TEST_NOTHROW( matrix.getLocalRowView(r,inds,vals) );
        TEST_COMPARE_ARRAYS( inds, tuple<LO>(r) );
        TEST_COMPARE_ARRAYS( vals, tuple<Scalar>(static_cast<Scalar>(3.0)) );
      }
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, CopiesAndViews, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // test that an exception is thrown when we exceed statically allocated memory
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = size(*comm);
    const size_t myImageID = rank(*comm);
    if (numImages < 2) return;
    // create a Map, one row per processor
    const size_t numLocal = 1;
    RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    GO myrowind = rmap->getGlobalElement(0);
    // specify the column map to control ordering
    // construct tridiagonal graph
    Array<GO> ginds;
    Array<LO> linds;
    if (myImageID==0) {
      DYN_ARRAYVIEW_TO_ARRAY( GO, ginds, tuple<GO>(myrowind,myrowind+1) );
      DYN_ARRAYVIEW_TO_ARRAY( LO, linds, tuple<LO>(0,1) );
    }
    else if (myImageID==numImages-1) {
      DYN_ARRAYVIEW_TO_ARRAY( GO, ginds , tuple<GO>(myrowind-1,myrowind) );
      DYN_ARRAYVIEW_TO_ARRAY( LO, linds , tuple<LO>(0,1) );
    }
    else {
      DYN_ARRAYVIEW_TO_ARRAY( GO, ginds , tuple<GO>(myrowind-1,myrowind,myrowind+1) );
      DYN_ARRAYVIEW_TO_ARRAY( LO, linds , tuple<LO>(0,1,2) );
    }
    Array<Scalar> vals(ginds.size(),ST::one());
    RCP<Map<LO,GO,Node> > cmap = rcp( new Map<LO,GO,Node>(INVALID,ginds(),0,comm,node) );
    for (int T=0; T<4; ++T) {
      ProfileType pftype = ( (T & 1) == 1 ) ? StaticProfile : DynamicProfile;
      OptimizeOption os  = ( (T & 2) == 2 ) ? DoOptimizeStorage : DoNotOptimizeStorage;
      MAT matrix(rmap,cmap, ginds.size(), pftype);   // only allocate as much room as necessary
      RowMatrix<Scalar,LO,GO,Node> &rowmatrix = matrix;
      Array<GO> GCopy(4); Array<LO> LCopy(4); Array<Scalar> SCopy(4);
      ArrayView<const GO> CGView; ArrayView<const LO> CLView; ArrayView<const Scalar> CSView;
      size_t numentries;
      // at this point, the graph has not allocated data as global or local, so we can do views/copies for either local or global
      matrix.getLocalRowCopy(0,LCopy,SCopy,numentries);
      matrix.getLocalRowView(0,CLView,CSView);
      matrix.getGlobalRowCopy(myrowind,GCopy,SCopy,numentries);
      matrix.getGlobalRowView(myrowind,CGView,CSView);
      // use multiple inserts: this illustrated an overwrite bug for column-map-specified graphs
      for (size_t j=0; j<(size_t)ginds.size(); ++j) {
        matrix.insertGlobalValues(myrowind,ginds(j,1),tuple(ST::one()));
      }
      TEST_EQUALITY( matrix.getNumEntriesInLocalRow(0), matrix.getCrsGraph()->getNumAllocatedEntriesInLocalRow(0) ); // test that we only allocated as much room as necessary
      // if static graph, insert one additional entry on my row and verify that an exception is thrown
      if (pftype == StaticProfile) {
        TEST_THROW( matrix.insertGlobalValues(myrowind,arrayView(&myrowind,1),tuple(ST::one())), std::runtime_error );
      }
      matrix.fillComplete(os);
      // check for throws and no-throws/values
      TEST_THROW( matrix.getGlobalRowView(myrowind,CGView,CSView), std::runtime_error );
      TEST_THROW( matrix.getLocalRowCopy(    0       ,LCopy(0,1),SCopy(0,1),numentries), std::runtime_error );
      TEST_THROW( matrix.getGlobalRowCopy(myrowind,GCopy(0,1),SCopy(0,1),numentries), std::runtime_error );
      //
      TEST_NOTHROW( matrix.getLocalRowView(0,CLView,CSView) );
      TEST_COMPARE_ARRAYS( CLView, linds );
      TEST_COMPARE_ARRAYS( CSView, vals  );
      //
      TEST_NOTHROW( matrix.getLocalRowCopy(0,LCopy,SCopy,numentries) );
      TEST_COMPARE_ARRAYS( LCopy(0,numentries), linds );
      TEST_COMPARE_ARRAYS( SCopy(0,numentries), vals  );
      //
      TEST_NOTHROW( matrix.getGlobalRowCopy(myrowind,GCopy,SCopy,numentries) );
      TEST_COMPARE_ARRAYS( GCopy(0,numentries), ginds );
      TEST_COMPARE_ARRAYS( SCopy(0,numentries), vals  );
      // 
      STD_TESTS(rowmatrix);
    }
    // All procs fail if any node fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TheEyeOfTruth, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();
    // create the identity matrix
    GO base = numLocal*myImageID;
    RCP<RowMatrix<Scalar,LO,GO,Node> > eye;
    {
      RCP<MAT> eye_crs = rcp(new MAT(map,1));
      for (size_t i=0; i<numLocal; ++i) {
        eye_crs->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<Scalar>(ST::one()));
      }
      TEST_EQUALITY_CONST( eye_crs->getProfileType() == DynamicProfile, true );
      eye_crs->fillComplete();
      eye = eye_crs;
    }
    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumRows()          , numLocal);
    TEST_EQUALITY(eye->getNodeNumCols()          , numLocal);
    TEST_EQUALITY(eye->getGlobalNumDiags() , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumDiags()     , numLocal);
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getNodeMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getRangeMap()) , true);
    // test the action
    mvres.randomize();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NonSquare, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t M = 3;
    const size_t P = 5;
    const size_t N = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create Maps
    // matrix is M*N-by-P
    //                  col
    //            0        1                  P-1
    //    0  [0        MN              ... (P-1)MN     ]
    //    .  [...      ...                 ...         ]
    //    0  [M-1      MN+M-1              (P-1)MN+M-1 ]
    //p   1  [M        MN+M                            ]
    //r   .  [...      ...                             ] = [A_ij], where A_ij = i+jMN
    //o   1  [2M-1     MN+2M-1                         ]
    //c   .  [...                                      ]
    //   N-1 [(N-1)M   MN+(N-1)(M-1)                   ]    
    //    .  [...      ...                             ]
    //   N-1 [MN-1     MN+MN-1                         ]
    // 
    // row map, range map is [0,M-1] [M,2M-1] [2M,3M-1] ... [MN-M,MN-1]
    // domain map will be map for X (lclmap)
    // 
    // input multivector X is not distributed:
    //               
    //   X = [  0    P    ...  (numVecs-1)P ]
    //       [ ...  ....  ...       ...     ] = [X_ji], where X_ij = i+jP
    //       [ P-1  2P-1  ...   numVecs*P-1 ]
    //
    // the result of the non-transpose multiplication should be 
    //                              P-1
    // (A*X)_ij = sum_k A_ik X_kj = sum (i+kMN)(k+jP) = jiP^2 + (i+jMNP)(P^2-P)/2 + MNP(P-1)(2P-1)/6
    //                              k=0
    // 
    // 
    // 
    const size_t numVecs  = 3;
    RCP<const Map<LO,GO,Node> > rowmap = createContigMapWithNode<LO,GO>(INVALID,M,comm,node);
    RCP<const Map<LO,GO,Node> > lclmap = createLocalMapWithNode<LO,GO,Node>(P,comm,node);

    // create the matrix
    MAT A(rowmap,P,DynamicProfile);
    for (GO i=0; i<static_cast<GO>(M); ++i) {
      for (GO j=0; j<static_cast<GO>(P); ++j) {
        A.insertGlobalValues( M*myImageID+i, tuple<GO>(j), tuple<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    TEST_EQUALITY_CONST( A.getProfileType() == DynamicProfile, true );
    A.fillComplete(lclmap,rowmap);
    // build the input multivector X
    MV X(lclmap,numVecs);
    for (GO i=0; i<static_cast<GO>(P); ++i) {
      for (GO j=0; j<static_cast<GO>(numVecs); ++j) {
        X.replaceGlobalValue(i,j,static_cast<Scalar>(i+j*P));
      }
    }
    // build the expected output multivector B
    MV Bexp(rowmap,numVecs), Bout(rowmap,numVecs);
    for (GO i=static_cast<GO>(myImageID*M); i<static_cast<GO>(myImageID*M+M); ++i) {
      for (GO j=0; j<static_cast<GO>(numVecs); ++j) {
        Bexp.replaceGlobalValue(i,j,static_cast<Scalar>(j*i*P*P + (i+j*M*N*P)*(P*P-P)/2 + M*N*P*(P-1)*(2*P-1)/6));
      }
    }
    // test the action
    Bout.randomize();
    A.apply(X,Bout);
    Bout.update(-ST::one(),Bexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    Bout.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Transpose, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // this is the same matrix as in test NonSquare, but we will apply the transpose
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t M = 3;
    const size_t P = 5;
    const size_t N = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create Maps
    // matrix is M*N-by-P
    //                  col
    //            0        1                  P-1
    //    0  [0        MN              ... (P-1)MN     ]
    //    .  [...      ...                 ...         ]
    //    0  [M-1      MN+M-1              (P-1)MN+M-1 ]
    //p   1  [M        MN+M                            ]
    //r   .  [...      ...                             ] = [A_ij], where A_ij = i+jMN
    //o   1  [2M-1     MN+2M-1                         ]
    //c   .  [...                                      ]
    //   N-1 [(N-1)M   MN+(N-1)(M-1)                   ]    
    //    .  [...      ...                             ]
    //   N-1 [MN-1     MN+MN-1                         ]
    // 
    // row map, range map is [0,M-1] [M,2M-1] [2M,3M-1] ... [MN-M,MN-1]
    // domain map will be a non-distributed map for a vector of length P
    // 
    // input multivector is 
    //  col      0            1         ...        numVecs-1
    //     0 [0         MN                    (numVecs-1)MN     ]
    // p   . [...       ...                   ...               ]
    // r   0 [M-1       MN+M-1                (numVecs-1)MN+M-1 ]
    // o   1 [M         MN+M                                    ]
    // c   . [...       ...                                     ] = [X_ij], where X_ij = i+jMN
    //     1 [2M-1      MN+2M-1                                 ]
    //     . [...       ...                                     ]
    //    N-1[(N-1)M    MN+(N-1)(M-1)                           ]
    //     . [...       ...                                     ]
    //    N-1[MN-1      MN+MN-1                                 ]
    //
    // output multivector is not-distributed
    // the result of the transpose multiplication should be 
    //              MN-1              MN-1
    // (A^T*X)_ij = sum_k A_ki X_kj = sum (k+iMN)(k+jMN) 
    //              k=0               k=0
    //   MN-1
    // = sum k(i+j)MN + ij(MN)(MN) + k^2 = (i+j)(MN)^2(MN-1)/2 + ij(MN)^3 + (MN)(MN-1)(2MN-1)/6
    //   k=0
    // 
    const size_t numVecs  = 3;
    RCP<const Map<LO,GO,Node> > rowmap = createContigMapWithNode<LO,GO>(INVALID,M,comm,node);
    RCP<const Map<LO,GO,Node> > lclmap = createLocalMapWithNode<LO,GO,Node>(P,comm,node);
    // create the matrix
    MAT A(rowmap,P);
    for (size_t i=0; i<M; ++i) {
      for (size_t j=0; j<P; ++j) {
        A.insertGlobalValues( static_cast<GO>(M*myImageID+i), tuple<GO>(j), tuple<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    TEST_EQUALITY_CONST( A.getProfileType() == DynamicProfile, true );
    A.fillComplete(lclmap,rowmap);
    out << "A: " << endl << A << endl;
    // build the input multivector X
    MV X(rowmap,numVecs);
    for (size_t i=myImageID*M; i<myImageID*M+M; ++i) {
      for (size_t j=0; j<numVecs; ++j) {
        X.replaceGlobalValue(i,j,static_cast<Scalar>( i + j*M*N ) );
      }
    }
    // build the expected output multivector B
    MV Bexp(lclmap,numVecs), Bout(lclmap,numVecs);
    for (size_t i=0; i<P; ++i) {
      for (size_t j=0; j<numVecs; ++j) {
        Bexp.replaceGlobalValue(i,j,static_cast<Scalar>( (i+j)*(M*N)*(M*N)*(M*N-1)/2 + i*j*(M*N)*(M*N)*(M*N) + (M*N)*(M*N-1)*(2*M*N-1)/6 ));
      }
    }
    // test the action
    Bout.randomize();
    A.apply(X,Bout,CONJ_TRANS);
    Bout.update(-ST::one(),Bexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    Bout.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, DomainRange, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages == 1) return;
    // create Maps
    // matrix is:
    //  proc
    //    0  [1 1               ]
    //    0  [1 1 1             ]
    //    1  [  1 1 1           ]
    //    1  [    1 1           ]
    //         ................
    //   n-1 [             1 1 1]
    //   n-1 [               1 1]
    //                              proc
    // input multivector will be      0  [0    2N   4N   6N    8N  ]
    //                                0  [1    2N+1 4N+1 6N+1  8N+1]
    //                                1  [2    2N+2 4N+2 6N+2  8N+2]
    //                                1  [3    2N+3 4N+3 6N+3  8N+3]                  (i,j) = 2*N*j+i
    //                                   [                         ]
    //                               n-1 [2N-2 4N-2 6N-2 8N-2 10N-2]
    //                               n-1 [2N-1 4N-1 6N-1 8N-1 10N-1]
    //                
    //                              proc                                                 
    // output multivector will be     0  [1    4N+1  8N+1 12N+1 16N+1]                    i=0:         (i,j)+(i+1,j) = 2Nj+i+2Nj+i+1 = 4Nj+2i+1 = 4Nj+1
    //                                1  [3    6N+3 12N+3 18N+3 24N+3]                 i=2n-1: (i-1,j)+(i,j)         = 2Nj+i-1+2Nj+i = 4Nj+2i-1 = 4Nj+2(2N-1)-1 = 4N(j+1)-3
    //                                2  [6    6N+6 12N+6 18N+6 24N+6]                   else: (i-1,j)+(i,j)+(i+1,j) = 2Nj+i-1+2Nj+i+2Nj+i+1 = 6Nj+3i
    //                                   [                           ]
    //                               n-1 [                           ]
    //                                0  [                           ]
    //                                1  [                           ]
    //                                   [                           ]
    //                               n-1 [4N-3 8N-3 12N-3 16N-3 20N-3]
    //
    // row map is [0,1]   [2,3]     [4,5]     etc
    // col map is [0,1,2] [1,2,3,4] [3,4,5,6] etc     (assembled by CrsMatrix, we construct one only for comparison)
    // domain map will be equal to the row map
    // range  map will be [0,np] [1,np+1] [2,np+2]
    const size_t numVecs  = 5;
    RCP<Map<LO,GO,Node> > rowmap = rcp( new Map<LO,GO,Node>(INVALID,tuple<GO>(2*myImageID,2*myImageID+1),0,comm,node) );
    RCP<Map<LO,GO,Node> > rngmap = rcp( new Map<LO,GO,Node>(INVALID,tuple<GO>(myImageID,numImages+myImageID),0,comm,node) );
    RCP<RowMatrix<Scalar,LO,GO,Node> > tri;
    {
      RCP<MAT> tri_crs = rcp(new MAT(rowmap,3) );
      Array<Scalar>  vals(3,ST::one());
      if (myImageID == 0) {
        Array<GO> cols( tuple<GO>(2*myImageID,2*myImageID+1,2*myImageID+2) );
        tri_crs->insertGlobalValues(2*myImageID  ,cols(0,2),vals(0,2));
        tri_crs->insertGlobalValues(2*myImageID+1,cols(0,3),vals(0,3));
      }
      else if (myImageID == numImages-1) {        
        Array<GO> cols( tuple<GO>(2*myImageID-1,2*myImageID,2*myImageID+1) );
        tri_crs->insertGlobalValues(2*myImageID  ,cols(0,3),vals(0,3));
        tri_crs->insertGlobalValues(2*myImageID+1,cols(1,2),vals(1,2));
      }
      else {
        Array<GO> cols( tuple<GO>(2*myImageID-1,2*myImageID,2*myImageID+1,2*myImageID+2) );
        tri_crs->insertGlobalValues(2*myImageID  ,cols(0,3),vals(0,3));
        tri_crs->insertGlobalValues(2*myImageID+1,cols(1,3),vals(0,3));
      }
      // call fillComplete(), specifying domain and range maps and requiring custom importer and exporter
      tri_crs->fillComplete(rowmap,rngmap);
      tri = tri_crs;
    }
    // test the properties
    TEST_EQUALITY(tri->getGlobalNumEntries()  , static_cast<size_t>(6*numImages-2));          
    TEST_EQUALITY(tri->getNodeNumEntries()      , (myImageID > 0 && myImageID < numImages-1) ? 6 : 5);
    TEST_EQUALITY(tri->getGlobalNumRows()      , static_cast<size_t>(2*numImages));
    TEST_EQUALITY(tri->getNodeNumRows()          , 2);
    TEST_EQUALITY(tri->getNodeNumCols()          , (myImageID > 0 && myImageID < numImages-1) ? 4 : 3);
    TEST_EQUALITY(tri->getGlobalNumDiags() , static_cast<size_t>(2*numImages));
    TEST_EQUALITY(tri->getNodeNumDiags()     , 2);
    TEST_EQUALITY(tri->getGlobalMaxNumRowEntries(), 3);
    TEST_EQUALITY(tri->getNodeMaxNumRowEntries()    , 3);
    TEST_EQUALITY(tri->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(tri->getRowMap()->isSameAs(*rowmap), true);
    TEST_EQUALITY_CONST(tri->getRangeMap()->isSameAs(*rngmap), true);
    TEST_EQUALITY_CONST(tri->getDomainMap()->isSameAs(*rowmap), true);
    // build the input and corresponding output multivectors
    MV mvin(rowmap,numVecs), mvout(rngmap,numVecs), mvexp(rngmap,numVecs);
    for (size_t j=0; j<numVecs; ++j) {
      mvin.replaceLocalValue(0,j,static_cast<Scalar>(j*2*numImages+2*myImageID  )); // entry (2*myImageID  ,j)
      mvin.replaceLocalValue(1,j,static_cast<Scalar>(j*2*numImages+2*myImageID+1)); // entry (2*myImageID+1,j)
      // entry (myImageID,j)
      if (myImageID==0) {
        mvexp.replaceLocalValue(0,j,static_cast<Scalar>(4*numImages*j+1));
      }                                                                    
      else {                                                               
        mvexp.replaceLocalValue(0,j,static_cast<Scalar>(6*numImages*j+3*myImageID));
      }                                                                    
      // entry (numImages+myImageID,j)
      if (myImageID==numImages-1) {                                        
        mvexp.replaceLocalValue(1,j,static_cast<Scalar>(4*numImages*(j+1)-3));
      }
      else {
        mvexp.replaceLocalValue(1,j,static_cast<Scalar>(6*numImages*j+3*(numImages+myImageID)));
      }
    }
    // test the action
    mvout.randomize();
    tri->apply(mvin,mvout);
    mvout.update(-ST::one(),mvexp,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvout.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TheEyeOfTruthDistAlloc, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();
    // create the identity matrix
    RCP<RowMatrix<Scalar,LO,GO,Node> > eye;
    {
      RCP<MAT> eye_crs = rcp(new MAT(map,1) );
      if (myImageID == 0) {
        for (int i=0; i<map->getGlobalNumEntries(); ++i) {
          eye_crs->insertGlobalValues(i,tuple<GO>(i),tuple<Scalar>(ST::one()));
        }
      }
      eye_crs->fillComplete();
      eye = eye_crs;
    }
    // test the properties
    TEST_EQUALITY(eye->getGlobalNumEntries()  , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumEntries()      , numLocal);
    TEST_EQUALITY(eye->getGlobalNumRows()      , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumRows()          , numLocal);
    TEST_EQUALITY(eye->getNodeNumCols()          , numLocal);
    TEST_EQUALITY(eye->getGlobalNumDiags() , numImages*numLocal);
    TEST_EQUALITY(eye->getNodeNumDiags()     , numLocal);
    TEST_EQUALITY(eye->getGlobalMaxNumRowEntries(), 1);
    TEST_EQUALITY(eye->getNodeMaxNumRowEntries()    , 1);
    TEST_EQUALITY(eye->getIndexBase()          , 0);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getColMap())   , true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getDomainMap()), true);
    TEST_EQUALITY_CONST(eye->getRowMap()->isSameAs(*eye->getRangeMap()) , true);
    // test the action
    mvres.randomize();
    eye->apply(mvrand,mvres);
    mvres.update(-ST::one(),mvrand,ST::one());
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, SimpleEigTest, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t ONE = OrdinalTraits<size_t>::one();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 2) return;
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,ONE,comm,node);
    // create a multivector ones(n,1)
    MV ones(map,ONE,false), threes(map,ONE,false);
    ones.putScalar(ST::one());
    /* create the following matrix:
       [2 1           ]
       [1 1 1         ]
       [  1 1 1       ]
       [   . . .      ]
       [     . . .    ]
       [       . . .  ]
       [         1 1 1]
       [           1 2]
     this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
    */
    size_t myNNZ;
    MAT A(map,3);
    if (myImageID == 0) {
      myNNZ = 2;
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*ST::one(), ST::one()));
      Array<GO> cols(tuple<GO>(myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      Array<Scalar> vals(tuple<Scalar>(ST::one(), static_cast<Scalar>(2)*ST::one()));
      Array<GO> cols(tuple<GO>(myImageID-1,myImageID));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    else {
      myNNZ = 3;
      Array<Scalar> vals(3,ST::one());
      Array<GO> cols(tuple<GO>(myImageID-1, myImageID, myImageID+1));
      A.insertGlobalValues(myImageID,cols(),vals());
    }
    A.fillComplete();
    // test the properties
    TEST_EQUALITY(A.getGlobalNumEntries()   , static_cast<size_t>(3*numImages-2));
    TEST_EQUALITY(A.getNodeNumEntries()       , myNNZ);
    TEST_EQUALITY(A.getGlobalNumRows()       , static_cast<size_t>(numImages));
    TEST_EQUALITY_CONST(A.getNodeNumRows()     , ONE);
    TEST_EQUALITY(A.getNodeNumCols()           , myNNZ);
    TEST_EQUALITY(A.getGlobalNumDiags()  , static_cast<size_t>(numImages));
    TEST_EQUALITY_CONST(A.getNodeNumDiags(), ONE);
    TEST_EQUALITY(A.getGlobalMaxNumRowEntries() , (numImages > 2 ? 3 : 2));
    TEST_EQUALITY(A.getNodeMaxNumRowEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getDomainMap()), true);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getRangeMap()) , true);
    // test the action
    threes.randomize();
    A.apply(ones,threes);
    // now, threes should be 3*ones
    threes.update(static_cast<Scalar>(-3)*ST::one(),ones,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    threes.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EmptyTriSolve, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Operator<Scalar,LO,GO,Node>  OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t numLocal = 13, numVecs = 7;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);

    /* Create a triangular matrix with no entries, for testing implicit diagonals.
      We test with Transpose and Non-Transpose application solve (these should be equivalent for the identity matrix)
    */
    
    MV X(map,numVecs), B(map,numVecs), Xhat(map,numVecs);
    X.setObjectLabel("X");
    B.setObjectLabel("B");
    Xhat.setObjectLabel("Xhat");
    X.randomize();
    for (size_t tnum=0; tnum < 2; ++tnum) {
      ETransp trans     = ((tnum & 1) == 1 ? CONJ_TRANS        : NO_TRANS);
      RCP<OP> ZeroIOp;
      {
        RCP<MAT> ZeroMat;
        // must explicitly provide the column map for implicit diagonals
        ZeroMat = rcp(new MAT(map,map,0));
        ZeroMat->fillComplete(DoOptimizeStorage);
        TEST_EQUALITY_CONST(ZeroMat->isUpperTriangular(), true);
        TEST_EQUALITY_CONST(ZeroMat->isLowerTriangular(), true);
        TEST_EQUALITY_CONST(ZeroMat->getGlobalNumDiags(), 0);
        ZeroIOp = createCrsMatrixSolveOp<Scalar>(ZeroMat.getConst());
      }
      X = B;
      Xhat.randomize();
      ZeroIOp->apply(B,Xhat,trans);
      //
      Xhat.update(-ST::one(),X,ST::one());
      Array<Mag> errnrms(numVecs), normsB(numVecs), zeros(numVecs, MT::zero());
      Xhat.norm2(errnrms());
      B.norm2(normsB());
      Mag maxBnrm = *std::max_element( normsB.begin(), normsB.end() );
      TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TriSolve, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Operator<Scalar,LO,GO,Node>  OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t numLocal = 13, numVecs = 7;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    Scalar SONE = static_cast<Scalar>(1.0);

    /* Create one of the following locally triangular matries:
     
    0  [1 2       ] 
    1  [  1 3     ] 
    .  [    .  .  ] = U
   n-2 [       1 n]
   n-1 [         1]

    0  [1           ] 
    1  [2 1         ] 
    .  [   .  .     ] = L
   n-2 [     n-1 1  ]
   n-1 [         n 1]

      Global matrices are diag(U,U,...,U) and diag(L,L,...,L)
    
      For each of these, we test with explicit and implicit unit diagonal, Transpose and Non-Transpose application, 1D or 2D storage
      Ultimately, that is 16 combinations:
      (Upper vs. Lower)  x  (Explicit vs. Implicit diagonal)  x  (Transpose vs. Non-Transpose)  x  (Optimized vs. Non-Optimzied storage)
    */
    
    MV X(map,numVecs), B(map,numVecs), Xhat(map,numVecs);
    X.setObjectLabel("X");
    B.setObjectLabel("B");
    Xhat.setObjectLabel("Xhat");
    X.randomize();
    for (size_t tnum=0; tnum < 16; ++tnum) {
      EUplo   uplo      = ((tnum & 1) == 1 ? UPPER_TRI         : LOWER_TRI);
      EDiag   diag      = ((tnum & 2) == 2 ? UNIT_DIAG         : NON_UNIT_DIAG);
      OptimizeOption os = ((tnum & 4) == 4 ? DoOptimizeStorage : DoNotOptimizeStorage);
      ETransp trans     = ((tnum & 8) == 8 ? CONJ_TRANS        : NO_TRANS);
      RCP<OP> AIOp;
      {
        RCP<MAT> AMat;
        if (diag == UNIT_DIAG) {
          // must explicitly specify the column map
          AMat = rcp(new MAT(map,map,2));
        }
        else {
          // can let the matrix compute a column map
          AMat = rcp(new MAT(map,2));
        }
        // fill the matrix
        if (uplo == UPPER_TRI) {
          if (diag == UNIT_DIAG) {
            for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
              if (gid == map->getMaxGlobalIndex()) {
                // do nothing
              }
              else {
                AMat->insertGlobalValues( gid, tuple<GO>(gid+1), tuple<Scalar>(static_cast<GO>(gid+2)) );
              }
            }
          }
          else {
            for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
              if (gid == map->getMaxGlobalIndex()) {
                AMat->insertGlobalValues( gid, tuple<GO>(gid), tuple<Scalar>(SONE) );
              }
              else {
                AMat->insertGlobalValues( gid, tuple<GO>(gid,gid+1), tuple<Scalar>(SONE,static_cast<GO>(gid+2)) );
              }
            }
          }
        }
        else { // uplo == LOWER_TRI
          if (diag == UNIT_DIAG) {
            for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
              if (gid == map->getMinGlobalIndex()) {
                // do nothing
              }
              else {
                AMat->insertGlobalValues( gid, tuple<GO>(gid-1), tuple<Scalar>(static_cast<GO>(gid+1)) );
              }
            }
          }
          else {
            for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
              if (gid == map->getMinGlobalIndex()) {
                AMat->insertGlobalValues( gid, tuple<GO>(gid), tuple<Scalar>(SONE) );
              }
              else {
                AMat->insertGlobalValues( gid, tuple<GO>(gid-1,gid), tuple<Scalar>(static_cast<GO>(gid+1),SONE) );
              }
            }
          }
        }
        AMat->fillComplete(os);
        TEST_EQUALITY(AMat->isUpperTriangular(), uplo == UPPER_TRI);
        TEST_EQUALITY(AMat->isLowerTriangular(), uplo == LOWER_TRI);
        TEST_EQUALITY(AMat->getGlobalNumDiags() == 0, diag == UNIT_DIAG);
        AIOp = createCrsMatrixSolveOp<Scalar>(AMat.getConst());
      }
      B.randomize();
      AIOp->apply(X,B,trans);
      if (diag == UNIT_DIAG) {
        // we want (I+A)*X -> B
        // A*X -> B needs to be augmented with X
        B.update(ST::one(),X,ST::one());
      }
      Xhat.randomize();
      AIOp->apply(B,Xhat,trans);
      //
      Xhat.update(-ST::one(),X,ST::one());
      Array<Mag> errnrms(numVecs), normsB(numVecs), zeros(numVecs, MT::zero());
      Xhat.norm2(errnrms());
      B.norm2(normsB());
      Mag maxBnrm = *std::max_element( normsB.begin(), normsB.end() );
      TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, FullMatrixTriDiag, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // do a FEM-type communication, then apply to a MultiVector containing the identity
    // this will check non-trivial communication and test multivector apply
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t ONE = OrdinalTraits<size_t>::one();
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numImages = comm->getSize();
    const size_t myImageID = comm->getRank();
    if (numImages < 3) return;
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,ONE,comm,node);

    // RCP<FancyOStream> fos = Teuchos::fancyOStream(rcp(&std::cout,false));
    
    /* Create the following matrix:
    0  [2 1       ]   [2 1]
    1  [1 4 1     ]   [1 2] + [2 1]
    2  [  1 4 1   ]           [1 2] + 
    3  [    1     ] = 
       [       4 1]
   n-1 [       1 2]
    */
    size_t myNNZ;
    MAT A(map,4);
    A.setObjectLabel("The Matrix");
    MV mveye(map,numImages), mvans(map,numImages), mvres(map,numImages,true);
    mveye.setObjectLabel("mveye");
    mvans.setObjectLabel("mvans");
    mvres.setObjectLabel("mvres");
    if (myImageID != numImages-1) { // last image assigns none
      Array<Scalar> vals(tuple<Scalar>(static_cast<Scalar>(2)*ST::one(),ST::one(),static_cast<Scalar>(2)*ST::one()));
      Array<GO> cols(tuple<GO>(myImageID,myImageID + 1));
      A.insertGlobalValues(myImageID  ,cols(),vals(0,2)); // insert [2 1]
      A.insertGlobalValues(myImageID+1,cols(),vals(1,2)); // insert [1 2]
    }
    // put one on the diagonal, setting mveye to the identity
    mveye.replaceLocalValue(0,myImageID,ST::one());
    // divine myNNZ and build multivector with matrix
    if (myImageID == 0) {
      myNNZ = 2;
      mvans.replaceLocalValue(0,0,static_cast<Scalar>(2));
      mvans.replaceLocalValue(0,1,static_cast<Scalar>(1));
    }
    else if (myImageID == numImages-1) {
      myNNZ = 2;
      mvans.replaceLocalValue(0,numImages-2,static_cast<Scalar>(1));
      mvans.replaceLocalValue(0,numImages-1,static_cast<Scalar>(2));
    }
    else {
      myNNZ = 3;
      mvans.replaceLocalValue(0,myImageID-1,static_cast<Scalar>(1));
      mvans.replaceLocalValue(0,myImageID  ,static_cast<Scalar>(4));
      mvans.replaceLocalValue(0,myImageID+1,static_cast<Scalar>(1));
    }
    A.fillComplete();

    // test the properties
    TEST_EQUALITY(A.getGlobalNumEntries()     , static_cast<size_t>(3*numImages-2));
    TEST_EQUALITY(A.getNodeNumEntries()       , myNNZ);
    TEST_EQUALITY(A.getGlobalNumRows()       , static_cast<size_t>(numImages));
    TEST_EQUALITY_CONST(A.getNodeNumRows()     , ONE);
    TEST_EQUALITY(A.getNodeNumCols()           , myNNZ);
    TEST_EQUALITY(A.getGlobalNumDiags()  , static_cast<size_t>(numImages));
    TEST_EQUALITY_CONST(A.getNodeNumDiags(), ONE);
    TEST_EQUALITY(A.getGlobalMaxNumRowEntries() , 3);
    TEST_EQUALITY(A.getNodeMaxNumRowEntries()     , myNNZ);
    TEST_EQUALITY_CONST(A.getIndexBase()     , 0);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getColMap())   , false);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getDomainMap()), true);
    TEST_EQUALITY_CONST(A.getRowMap()->isSameAs(*A.getRangeMap()) , true);
    // test the action
    A.apply(mveye,mvres);
    mvres.update(-ST::one(),mvans,ST::one());
    Array<Mag> norms(numImages), zeros(numImages,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, AlphaBetaMultiply, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef  Operator<Scalar,LO,GO,Node> OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t THREE = 3;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t myImageID = comm->getRank();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,THREE,comm,node);

    /* Create the identity matrix, three rows per proc */
    RCP<OP> AOp;
    {
      RCP<MAT> A = rcp(new MAT(map,1));
      A->insertGlobalValues(3*myImageID,  tuple<GO>(3*myImageID  ), tuple<Scalar>(ST::one()) );
      A->insertGlobalValues(3*myImageID+1,tuple<GO>(3*myImageID+1), tuple<Scalar>(ST::one()) );
      A->insertGlobalValues(3*myImageID+2,tuple<GO>(3*myImageID+2), tuple<Scalar>(ST::one()) );
      A->fillComplete(DoOptimizeStorage);
      AOp = A;
    }
    MV X(map,1), Y(map,1), Z(map,1);
    const Scalar alpha = ST::random(),
                  beta = ST::random();
    X.randomize();
    Y.randomize();
    // Z = alpha*X + beta*Y
    Z.update(alpha,X,beta,Y,ST::zero());
    // test the action: Y = alpha*I*X + beta*Y = alpha*X + beta*Y = Z
    AOp->apply(X,Y,NO_TRANS,alpha,beta);
    // Z -= Y  -> zero
    Z.update(-ST::one(),Y,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    Z.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, MultiplyOp, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef CrsMatrix<int,LO,GO,Node> IntMAT;
    typedef  Operator<Scalar,LO,GO,Node> OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t THREE = 3;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t myImageID = comm->getRank();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,THREE,comm,node);

    /* Create the integer identity matrix, three rows per proc, wrapped in a Op<Scalar>  */
    RCP<OP> AOp;
    {
      RCP<IntMAT> A = rcp(new IntMAT(map,1));
      A->insertGlobalValues(3*myImageID,  tuple<GO>(3*myImageID  ), tuple<int>(1) );
      A->insertGlobalValues(3*myImageID+1,tuple<GO>(3*myImageID+1), tuple<int>(1) );
      A->insertGlobalValues(3*myImageID+2,tuple<GO>(3*myImageID+2), tuple<int>(1) );
      A->fillComplete(DoOptimizeStorage);
      AOp = createCrsMatrixMultiplyOp<Scalar>(A.getConst());
    }
    MV X(map,1), Y(map,1), Z(map,1);
    X.randomize();
    Y.randomize();
    // Z = X + Y
    Z.update(ST::one(),X,ST::one(),Y,ST::zero());
    // test the action: Y = I*X + Y = X + Y == Z
    AOp->apply(X,Y,NO_TRANS,ST::one(),ST::one());
    // Z -= Y  -> zero
    Z.update(-ST::one(),Y,ST::one());
    Array<Mag> norms(1), zeros(1,MT::zero());
    Z.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, BadGID, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    // what happens when we call CrsMatrix::insertGlobalValues() for a row that isn't on the Map?
    typedef ScalarTraits<Scalar> ST;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t myImageID = comm->getRank();
    const size_t numImages = comm->getSize();
    // create a Map
    const size_t numLocal = 10;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    {
      // create the matrix
      MAT A(map,1);
      // add an entry off the map: row too high
      // this will only be off the map for the last node, for the others it will induce communication
      A.insertGlobalValues(map->getMaxGlobalIndex()+1,tuple<GO>(map->getIndexBase()),tuple<Scalar>(ST::one()));
      TEST_THROW(A.fillComplete(), std::runtime_error);
    }
    {
      // create the matrix
      MAT A(map,1);
      // add an entry off the map: row too high
      // this will only be off the map for the last node, for the others there is nothing
      if (myImageID == numImages-1) {
        A.insertGlobalValues(map->getMaxAllGlobalIndex()+1,tuple<GO>(map->getIndexBase()),tuple<Scalar>(ST::one()));
      }
      TEST_THROW(A.fillComplete(), std::runtime_error);
    }
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ZeroMatrix, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal = 10;
    const size_t numVecs  = 5;
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);
    MV mvrand(map,numVecs,false), mvres(map,numVecs,false);
    mvrand.randomize();
    // create the zero matrix
    MAT zero(map,0);
    zero.fillComplete();
    mvres.randomize();
    zero.apply(mvrand,mvres);
    Array<Mag> norms(numVecs), zeros(numVecs,MT::zero());
    mvres.norm2(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,MT::zero());
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Typedefs, LO, GO, Scalar, Node )
  {
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef typename MAT::scalar_type         scalar_type;
    typedef typename MAT::local_ordinal_type  local_ordinal_type;
    typedef typename MAT::global_ordinal_type global_ordinal_type;
    typedef typename MAT::node_type           node_type;
    typedef typename MAT::mat_vec_type        mat_vec_type;
    typedef typename MAT::mat_solve_type      mat_solve_type;
    TEST_EQUALITY_CONST( (is_same< scalar_type         , Scalar >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< local_ordinal_type  , LO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< global_ordinal_type , GO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< node_type           , Node   >::value) == true, true );
    typedef RowMatrix<Scalar,LO,GO,Node> RMAT;
    typedef typename RMAT::scalar_type         rmat_scalar_type;
    typedef typename RMAT::local_ordinal_type  rmat_local_ordinal_type;
    typedef typename RMAT::global_ordinal_type rmat_global_ordinal_type;
    typedef typename RMAT::node_type           rmat_node_type;
    TEST_EQUALITY_CONST( (is_same< rmat_scalar_type         , Scalar >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rmat_local_ordinal_type  , LO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rmat_global_ordinal_type , GO     >::value) == true, true );
    TEST_EQUALITY_CONST( (is_same< rmat_node_type           , Node   >::value) == true, true );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, ActiveFill, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    const GO INVALID = OrdinalTraits<GO>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,1,comm,node);
    const Scalar SZERO = ScalarTraits<Scalar>::zero();
    {
      MAT matrix(map,map,0,DynamicProfile);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      matrix.insertLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) );
      //
      matrix.fillComplete(DoNotOptimizeStorage);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      TEST_THROW( matrix.insertLocalValues ( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      TEST_THROW( matrix.replaceLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      TEST_THROW( matrix.sumIntoLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ), std::runtime_error );
      TEST_THROW( matrix.setAllToScalar(SZERO),                                   std::runtime_error );
      TEST_THROW( matrix.scale(SZERO),                                            std::runtime_error );  
      TEST_THROW( matrix.globalAssemble(),                                        std::runtime_error );
      TEST_THROW( matrix.fillComplete(),                                          std::runtime_error );
    }
    {
      MAT matrix(map,map,0,DynamicProfile);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      matrix.insertLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) );
      //
      matrix.fillComplete(DoNotOptimizeStorage);
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
      //
      matrix.resumeFill();
      TEST_EQUALITY_CONST( matrix.isFillActive(),   true );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), false );
      TEST_NOTHROW( matrix.insertLocalValues ( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.replaceLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.sumIntoLocalValues( 0, tuple<LO>(0), tuple<Scalar>(0) ) );
      TEST_NOTHROW( matrix.setAllToScalar(SZERO)                                   );
      TEST_NOTHROW( matrix.scale(SZERO)                                            );
      TEST_NOTHROW( matrix.globalAssemble()                                        );
      //
      TEST_NOTHROW( matrix.fillComplete()                        );
      TEST_EQUALITY_CONST( matrix.isFillActive(),   false );
      TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node > Scalar 
    getNorm(RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > matrix)
  {
    Scalar mySum = ScalarTraits<Scalar>::zero();
    Array<LocalOrdinal> inds(matrix->getNodeMaxNumRowEntries());
    Array<Scalar> vals(matrix->getNodeMaxNumRowEntries());
    for(LocalOrdinal i =0; ((size_t)i)<matrix->getNodeNumRows(); ++i){
      size_t numRowEnts = matrix->getNumEntriesInLocalRow(i);
      ArrayView<const LocalOrdinal> indsView = inds();
      ArrayView<const Scalar> valsView = vals();
      matrix->getLocalRowView(i, indsView, valsView);
      for(size_t j=0; ((size_t)j)<numRowEnts; ++j){
        mySum += valsView[j]*valsView[j];
      }
    }
    Scalar totalSum = 0;
    Teuchos::reduceAll(*(matrix->getComm()), Teuchos::REDUCE_SUM, 1, &mySum, &totalSum);
    return ScalarTraits<Scalar>::squareroot(totalSum);

  }

  //Constructs to Tridiagonal Matricies and scales them by a vector.
  //One on the left and one on the right. Then compares the result to the known
  //correct matrix
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, LeftRightScale, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef ScalarTraits<Scalar> ST;
    typedef OrdinalTraits<LO> LOT;
    typedef Vector<Scalar,LO,GO,Node> VEC;
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    int numProcs = comm->getSize();
    int myRank = comm->getRank();

    global_size_t numGlobal = 4*numProcs;
    RCP<const Map<LO,GO,Node> > map = createUniformContigMapWithNode<LO,GO>(numGlobal,comm,node);
    RCP<VEC> vector = createVector<Scalar, LO, GO, Node>(map);
    vector->putScalar(2);
   
    RCP<MAT> matrix = createCrsMatrix<Scalar, LO, GO, Node>(map, 3);
    RCP<MAT> matrix2= createCrsMatrix<Scalar, LO, GO, Node>(map, 3);
    RCP<MAT> answerMatrix = createCrsMatrix<Scalar, LO, GO, Node>(map, 3);

    Array<Scalar> vals = tuple<Scalar>(1,2,3);
    Array<Scalar> answerVals = tuple<Scalar>(2,4,6);
    Array<GO> cols(3,0);
    for(
      GO i = Teuchos::as<GO>(myRank)*4; 
      i<(Teuchos::as<GO>(myRank)*4)+4; 
      ++i)
    {
      if(i==0){
        cols = tuple<GO>(0,1);
        matrix->insertGlobalValues(i, cols(), vals(1,2));
        matrix2->insertGlobalValues(i, cols(), vals(1,2));
        answerMatrix->insertGlobalValues(i, cols(), answerVals(1,2));
      }
      else if(i==(Teuchos::as<GO>(numProcs-1)*4)+3){
        cols = tuple<GO>(numGlobal-2, numGlobal-1);
        matrix->insertGlobalValues(i, cols(), vals(0,2));
        matrix2->insertGlobalValues(i, cols(), vals(0,2));
        answerMatrix->insertGlobalValues(i, cols(), answerVals(0,2));
      }
      else{
        cols = tuple<GO>(i-1,i,i+1);
        matrix->insertGlobalValues(i, cols(), vals());
        matrix2->insertGlobalValues(i, cols(), vals());
        answerMatrix->insertGlobalValues(i, cols(), answerVals());
      }
    }

    matrix->fillComplete();
    matrix2->fillComplete();
    answerMatrix->fillComplete();
    matrix->leftScale(*vector);
    matrix2->rightScale(*vector);

    RCP<MAT> diffMat1 = createCrsMatrix<Scalar, LO, GO, Node>(map,3);
    RCP<MAT> diffMat2 = createCrsMatrix<Scalar, LO, GO, Node>(map,3);
    Scalar sOne = ScalarTraits<Scalar>::one();
    Tpetra::MatrixMatrix::Add(*matrix, false, -sOne, *answerMatrix, false, sOne, diffMat1);
    diffMat1->fillComplete();
    Tpetra::MatrixMatrix::Add(*matrix2, false, -sOne, *answerMatrix, false, sOne, diffMat2);
    diffMat2->fillComplete();
    Scalar epsilon1 = getNorm(diffMat1)/getNorm(answerMatrix);
    Scalar epsilon2 = getNorm(diffMat2)/getNorm(answerMatrix);
    TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon1), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon1), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon2), <, 1e-10)
    TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon2), <, 1e-10)
  }




// 
// INSTANTIATIONS
//


typedef std::complex<float>  ComplexFloat;
typedef std::complex<double> ComplexDouble;

#define UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TheEyeOfTruth, LO, GO, SCALAR, NODE )  \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ZeroMatrix   , LO, GO, SCALAR, NODE )  \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, BadCalls     , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, SimpleEigTest, LO, GO, SCALAR, NODE )  \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, BadGID       , LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, FullMatrixTriDiag, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, DomainRange, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NonSquare, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Transpose, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithGraph, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithColMap, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EmptyFillComplete, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, AdvancedGraphUsage, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, WithGraph_replaceLocal, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ExceedStaticAlloc, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, MultipleFillCompletes, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NonlocalAfterResume, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, CopiesAndViews, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TriSolve, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, AlphaBetaMultiply, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, MultiplyOp, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EmptyTriSolve, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, ActiveFill, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Typedefs,      LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, LeftRightScale,      LO, GO, SCALAR, NODE ) 


#define UNIT_TEST_SERIALNODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, SerialNode )

typedef Kokkos::DefaultNode::DefaultNodeType DefaultNode;
#define UNIT_TEST_DEFAULTNODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, DefaultNode )

#ifdef HAVE_KOKKOS_TBB
#define UNIT_TEST_TBBNODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, TBBNode )
#else
#define UNIT_TEST_TBBNODE(LO, GO, SCALAR)
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
#define UNIT_TEST_TPINODE(LO, GO, SCALAR) \
      UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, SCALAR, TPINode )
#else
#define UNIT_TEST_TPINODE(LO, GO, SCALAR)
#endif

// don't test Kokkos node for MPI builds, because we probably don't have multiple GPUs per node
#if defined(HAVE_KOKKOS_THRUST) && !defined(HAVE_TPETRA_MPI)
// float
#if defined(HAVE_KOKKOS_CUDA_FLOAT)
#  define UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, float, ThrustGPUNode )
#else
#  define UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO)
#endif
// double
#if defined(HAVE_KOKKOS_CUDA_DOUBLE)
#  define UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, double, ThrustGPUNode )
#else
#  define UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO)
#endif
// complex<float>
#if defined(HAVE_KOKKOS_CUDA_COMPLEX_FLOAT)
#  define UNIT_TEST_THRUSTGPUNODE_COMPLEX_FLOAT(LO, GO) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, ComplexFloat, ThrustGPUNode )
#else
#  define UNIT_TEST_THRUSTGPUNODE_COMPLEX_FLOAT(LO, GO)
#endif
// complex<double>
#if defined(HAVE_KOKKOS_CUDA_COMPLEX_DOUBLE)
#  define UNIT_TEST_THRUSTGPUNODE_COMPLEX_DOUBLE(LO, GO) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR_NODE( LO, GO, ComplexDouble, ThrustGPUNode )
#else
#  define UNIT_TEST_THRUSTGPUNODE_COMPLEX_DOUBLE(LO, GO)
#endif
#else
// none
# define UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO)
# define UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO)
# define UNIT_TEST_THRUSTGPUNODE_COMPLEX_FLOAT(LO, GO)
# define UNIT_TEST_THRUSTGPUNODE_COMPLEX_DOUBLE(LO, GO)
#endif

#define UNIT_TEST_ALLCPUNODES(LO, GO, SCALAR) \
    UNIT_TEST_DEFAULTNODE(LO, GO, SCALAR) 

#define UNIT_TEST_FLOAT(LO, GO) \
    UNIT_TEST_ALLCPUNODES(LO, GO, float) \
    UNIT_TEST_THRUSTGPUNODE_FLOAT(LO, GO)

#define UNIT_TEST_DOUBLE(LO, GO) \
    UNIT_TEST_ALLCPUNODES(LO, GO, double) \
    UNIT_TEST_THRUSTGPUNODE_DOUBLE(LO, GO)

#define UNIT_TEST_COMPLEX_FLOAT(LO, GO) \
    UNIT_TEST_ALLCPUNODES(LO, GO, ComplexFloat) \
    UNIT_TEST_THRUSTGPUNODE_COMPLEX_FLOAT(LO, GO)

#define UNIT_TEST_COMPLEX_DOUBLE(LO, GO) \
    UNIT_TEST_ALLCPUNODES(LO, GO, ComplexDouble) \
    UNIT_TEST_THRUSTGPUNODE_COMPLEX_DOUBLE(LO, GO)

#if defined(HAVE_TPETRA_INST_DOUBLE)
  UNIT_TEST_DOUBLE(int, int)
#endif

# if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
    UNIT_TEST_COMPLEX_FLOAT(int, int)
# endif 

}
