// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace {

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::tuple;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::CrsMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createVector;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

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
// UNIT TEST(S)
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NonlocalAfterResume, LO, GO, Scalar, Node )
{
  using std::cerr;
  using std::endl;

  // test that an exception is thrown when we exceed statically allocated memory
  typedef ScalarTraits<Scalar> ST;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  // get a comm
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  const size_t numImages = size(*comm);
  const size_t myImageID = rank(*comm);

  comm->barrier ();
  if (myImageID == 0) {
    std::ostringstream os;
    os << "=== Tpetra::CrsMatrix nonlocal-after-resume test ===" << endl;
    cerr << os.str ();
  }
  comm->barrier ();

  const GO numLocal = 5;
  {
    std::ostringstream os;
    os << "  Proc " << myImageID << ": Create row Map with " << numLocal
       << " rows per process" << endl;
    cerr << os.str ();
  }
  RCP<const Map<LO,GO,Node> > rmap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);

  {
    std::ostringstream os;
    os << "  Proc " << myImageID << ": "
       << "Create a column Map with super- and sub-diagonal blocks" << endl;
    cerr << os.str ();
  }
  RCP<const Map<LO,GO,Node> > cmap;
  {
    Array<GO> cols;
    for (GO c = rmap->getMinGlobalIndex (); c <= rmap->getMaxGlobalIndex (); ++c) {
      cols.push_back (c);
    }
    if (rmap->getMinGlobalIndex() >= rmap->getMinAllGlobalIndex() + numLocal) {
      for (GO c = rmap->getMinGlobalIndex()-numLocal; c < rmap->getMinGlobalIndex(); ++c) {
        cols.push_back(c);
      }
    }
    if (rmap->getMaxGlobalIndex()+numLocal <= rmap->getMaxAllGlobalIndex()) {
      for (GO c = rmap->getMaxGlobalIndex()+1; c <= rmap->getMaxGlobalIndex()+numLocal; ++c) {
        cols.push_back(c);
      }
    }
    cmap = createNonContigMapWithNode<LO,GO,Node>(cols(), comm);
  }

  comm->barrier ();
  if (myImageID == 0) {
    std::ostringstream os;
    os << "  GLOBAL: Created column Map" << endl;
    cerr << os.str ();
  }
  comm->barrier ();

  {
    {
      std::ostringstream os;
      os << "  Proc " << myImageID << ": Insert diagonal entries" << endl;
      cerr << os.str ();
    }
    //----------------------------------------------------------------------
    // put in diagonal, locally
    //----------------------------------------------------------------------
    Tpetra::CrsMatrix<Scalar,LO,GO,Node> matrix(rmap,cmap,3);
    for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
      matrix.insertGlobalValues(r,tuple(r),tuple(ST::one()));
    }

    {
      std::ostringstream os;
      os << "  Proc " << myImageID << ": Fill-complete the matrix" << endl;
      cerr << os.str ();
    }
    // fill, but do not pack, because we will add new entries below
    RCP<ParameterList> params = parameterList();
    params->set("Optimize Storage",false);
    params->set("compute global constants",true);
    TEST_NOTHROW       ( matrix.fillComplete( params ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(),      true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), numLocal*numImages );
    TEST_EQUALITY      ( matrix.getLocalNumEntries(),   (size_t)numLocal   );

    comm->barrier ();
    if (myImageID == 0) {
      std::ostringstream os;
      os << "  GLOBAL: Done with first fillComplete" << endl;
      cerr << os.str ();
    }
    comm->barrier ();

    {
      std::ostringstream os;
      os << "  Proc " << myImageID << ": Insert super-diagonal entries" << endl;
      cerr << os.str ();
    }
    //----------------------------------------------------------------------
    // add super-diagonal, non-locally
    //----------------------------------------------------------------------
    // because fillComplete() was called above, we must call resumeFill() before adding new entries
    matrix.resumeFill();
    if (rmap->getMinGlobalIndex()+numLocal < rmap->getMaxAllGlobalIndex()) {
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        matrix.insertGlobalValues(r+numLocal,tuple(r),tuple(ST::one()));
      }
    }
    {
      std::ostringstream os;
      os << "  Proc " << myImageID << ": Fill-complete the matrix" << endl;
      cerr << os.str ();
    }
    // fill, but do not pack, because we will add new entries below
    params->set("Optimize Storage",false);
    params->set("compute global constants",true);
    TEST_NOTHROW       ( matrix.fillComplete( params ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(),      true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), false );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), 2*numLocal*numImages-numLocal );
    {
      size_t expected = numLocal;
      if (myImageID > 0) expected += numLocal; // super-diagonal
      TEST_EQUALITY( matrix.getLocalNumEntries(), expected );
    }

    comm->barrier ();
    if (myImageID == 0) {
      std::ostringstream os;
      os << "  GLOBAL: Done with second fillComplete" << endl;
      cerr << os.str ();
    }
    comm->barrier ();

    {
      std::ostringstream os;
      os << "  Proc " << myImageID << ": Insert sub-diagonal entries" << endl;
      cerr << os.str ();
    }
    //----------------------------------------------------------------------
    // add sub-diagonal block, non-locally
    //----------------------------------------------------------------------
    // because fillComplete() was called above, we must call resumeFill() before adding new entries
    matrix.resumeFill();
    if (rmap->getMinGlobalIndex() >= rmap->getMinAllGlobalIndex()+numLocal) {
      for (GO r=rmap->getMinGlobalIndex(); r <= rmap->getMaxGlobalIndex(); ++r) {
        matrix.insertGlobalValues(r-numLocal,tuple(r),tuple(ST::one()));
      }
    }
    {
      std::ostringstream os;
      os << "  Proc " << myImageID << ": Fill-complete the matrix (with "
         << "optimized storage)" << endl;
      cerr << os.str ();
    }
    // fill; it is okay to pack now
    params->set("Optimize Storage",true);
    params->set("compute global constants",true);
    TEST_NOTHROW       ( matrix.fillComplete( params ) );
    TEST_EQUALITY_CONST( matrix.isFillComplete(), true );
    TEST_EQUALITY_CONST( matrix.isStorageOptimized(), true );
    TEST_EQUALITY      ( matrix.getGlobalNumEntries(), 3*numLocal*numImages-2*numLocal );
    {
      size_t expected = numLocal;
      if (myImageID > 0)           expected += numLocal; // super-diagonal
      if (myImageID < numImages-1) expected += numLocal; // sub-diagonal
      TEST_EQUALITY( matrix.getLocalNumEntries(), expected );
    }

    comm->barrier ();
    if (myImageID == 0) {
      std::ostringstream os;
      os << "  GLOBAL: Done with third fillComplete" << endl;
      cerr << os.str ();
    }
    comm->barrier ();
  }
  // All procs fail if any process fails
  int globalSuccess_int = -1;
  Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );

  if (myImageID == 0) {
    std::ostringstream os;
    os << "=== Done with test (globally) ===" << endl;
    cerr << os.str ();
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NonlocalAfterResume, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
}
