// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_ConfigDefs.hpp>

#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#endif
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Exceptions.hpp>

//#include <XpetraExt_MatrixMatrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>

//#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

#ifdef XPETRA_TEST_USE_LONGLONG_GO
#define MatrixMarketFileToCrsMatrix MatrixMarketFileToCrsMatrix64
#endif

namespace {

using std::find;
using std::sort;
using Teuchos::arcp;
using Teuchos::Array;
using Teuchos::as;
using Teuchos::broadcast;
using Teuchos::Comm;
using Teuchos::OrdinalTraits;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Tuple;
using Teuchos::tuple;

using Xpetra::CrsMatrix;
using Xpetra::DefaultPlatform;
using Xpetra::Map;
using Xpetra::Matrix;

using Xpetra::viewLabel_t;

bool testMpi         = true;
double errorTolSlack = 1e+1;

#if (defined(HAVE_XPETRA_EPETRA) && defined(HAVE_XPETRA_EPETRAEXT)) || (defined(HAVE_XPETRA_TPETRA))
RCP<const Comm<int> > getDefaultComm() {
  if (testMpi) {
    return DefaultPlatform::getDefaultPlatform().getComm();
  }
  return rcp(new Teuchos::SerialComm<int>());
}
#endif

/////////////////////////////////////////////////////

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results");
}

//
// UNIT TESTS
//

/// unit test for matrix-matrix multiplication (both for Epetra and Tpetra)
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, Multiply_Epetra, M, MA, Scalar, LO, GO, Node) {
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, Multiply_Epetra64, M, MA, Scalar, LO, GO, Node) {
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, Multiply_Tpetra, M, MA, Scalar, LO, GO, Node) {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  // The matrix reader does not work with complex scalars. Skip all tests then.
  return;
#endif
#ifdef HAVE_XPETRA_TPETRA
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::CrsMatrix<Scalar, LO, GO, Node> CrsMatrixClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrapClass;

  // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  // yAB->describe(*fos, Teuchos::VERB_EXTREME);

  {  // Tpetra test

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    Xpetra::UnderlyingLib lib  = Xpetra::UseTpetra;

    // define map
    LO nEle                       = 6;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    // read in matrices
    typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > reader_type;

    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > tpA    = reader_type::readSparseFile("A.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > tpB    = reader_type::readSparseFile("B.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > tpAB   = reader_type::readSparseFile("AB.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > tpAtB  = reader_type::readSparseFile("AtB.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > tpABt  = reader_type::readSparseFile("ABt.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> > tpAtBt = reader_type::readSparseFile("AtBt.mat", comm);

    // transform to Xpetra
    Teuchos::RCP<CrsMatrixClass> xAmat    = Teuchos::rcp(new MA(tpA));
    Teuchos::RCP<CrsMatrixClass> xBmat    = Teuchos::rcp(new MA(tpB));
    Teuchos::RCP<CrsMatrixClass> xABmat   = Teuchos::rcp(new MA(tpAB));
    Teuchos::RCP<CrsMatrixClass> xAtBmat  = Teuchos::rcp(new MA(tpAtB));
    Teuchos::RCP<CrsMatrixClass> xABtmat  = Teuchos::rcp(new MA(tpABt));
    Teuchos::RCP<CrsMatrixClass> xAtBtmat = Teuchos::rcp(new MA(tpAtBt));

    Teuchos::RCP<MatrixClass> xA    = Teuchos::rcp(new CrsMatrixWrapClass(xAmat));
    Teuchos::RCP<MatrixClass> xB    = Teuchos::rcp(new CrsMatrixWrapClass(xBmat));
    Teuchos::RCP<MatrixClass> xAB   = Teuchos::rcp(new CrsMatrixWrapClass(xABmat));
    Teuchos::RCP<MatrixClass> xAtB  = Teuchos::rcp(new CrsMatrixWrapClass(xAtBmat));
    Teuchos::RCP<MatrixClass> xABt  = Teuchos::rcp(new CrsMatrixWrapClass(xABtmat));
    Teuchos::RCP<MatrixClass> xAtBt = Teuchos::rcp(new CrsMatrixWrapClass(xAtBtmat));

    Teuchos::RCP<MatrixClass> yAB   = Teuchos::rcp(new CrsMatrixWrapClass(map, 10));
    Teuchos::RCP<MatrixClass> yAtB  = Teuchos::rcp(new CrsMatrixWrapClass(map, 10));
    Teuchos::RCP<MatrixClass> yABt  = Teuchos::rcp(new CrsMatrixWrapClass(map, 10));
    Teuchos::RCP<MatrixClass> yAtBt = Teuchos::rcp(new CrsMatrixWrapClass(map, 10));

    Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply(
        *xA,
        false,
        *xB,
        false,
        *yAB);
    TEUCHOS_TEST_EQUALITY(xAB->getFrobeniusNorm(), yAB->getFrobeniusNorm(), out, success);
    TEUCHOS_TEST_EQUALITY(xAB->getLocalNumEntries(), yAB->getLocalNumEntries(), out, success);

    Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply(
        *xA,
        true,
        *xB,
        false,
        *yAtB);
    TEUCHOS_TEST_EQUALITY(xAtB->getFrobeniusNorm(), yAtB->getFrobeniusNorm(), out, success);
    TEUCHOS_TEST_EQUALITY(xAtB->getLocalNumEntries(), yAtB->getLocalNumEntries(), out, success);

    Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply(
        *xA,
        false,
        *xB,
        true,
        *yABt);
    TEUCHOS_TEST_EQUALITY(xABt->getFrobeniusNorm(), yABt->getFrobeniusNorm(), out, success);
    TEUCHOS_TEST_EQUALITY(xABt->getLocalNumEntries(), yABt->getLocalNumEntries(), out, success);

    Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply(
        *xA,
        true,
        *xB,
        true,
        *yAtBt);
    TEUCHOS_TEST_EQUALITY(xAtBt->getFrobeniusNorm(), yAtBt->getFrobeniusNorm(), out, success);
    TEUCHOS_TEST_EQUALITY(xAtBt->getLocalNumEntries(), yAtBt->getLocalNumEntries(), out, success);
  }  // end Tpetra test
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, BlockCrs, M, MB, Scalar, LO, GO, Node) {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  // The matrix reader does not work with complex scalars. Skip all tests then.
  return;
#endif
#ifdef HAVE_XPETRA_TPETRA
  typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
  typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;

  typedef Xpetra::CrsMatrix<Scalar, LO, GO, Node> CrsMatrixClass;
  typedef Xpetra::TpetraBlockCrsMatrix<Scalar, LO, GO, Node> BlockCrsMatrixClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrapClass;
  using helpers = Xpetra::Helpers<Scalar, LO, GO, Node>;

  RCP<const Comm<int> > comm = getDefaultComm();
  const GO INVALID           = Teuchos::OrdinalTraits<GO>::invalid();

  const size_t numLocalMeshPoints = 12;
  const GO indexBase              = 1;
  // mfh 16 May 2014: Tpetra::CrsGraph still needs the row Map as an
  // RCP.  Later interface changes will let us pass in the Map by
  // const reference and assume view semantics.
  RCP<const map_type> meshRowMapPtr =
      rcp(new map_type(INVALID, numLocalMeshPoints, indexBase, comm));

  // Test w/ an empty graph
  const LO blockSize = 4;
  graph_type graph(meshRowMapPtr, 0);
  graph.fillComplete();

  // Make the matrix (Tpetra)
  RCP<BCM> blockMat            = rcp(new BCM(graph, blockSize));
  RCP<CrsMatrixClass> bmc      = rcp(new BlockCrsMatrixClass(blockMat));
  RCP<CrsMatrixWrapClass> wrap = rcp(new CrsMatrixWrapClass(bmc));
  RCP<MatrixClass> mat         = wrap;

  // Now for the checks
  TEUCHOS_TEST_EQUALITY(helpers::isTpetraBlockCrs(mat), true, out, success);
  TEUCHOS_TEST_EQUALITY(helpers::isTpetraCrs(mat), false, out, success);

#endif
}

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                                  \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N;              \
  typedef typename Xpetra::TpetraCrsMatrix<S, LO, GO, N> MA##S##LO##GO##N; \
  typedef typename Xpetra::TpetraBlockCrsMatrix<S, LO, GO, N> MB##S##LO##GO##N;

#endif

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S, LO, GO, N)                                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, Multiply_Tpetra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, BlockCrs, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, Multiply_Epetra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

// List of tests which run only with Epetra64
#define XP_EPETRA64_MATRIX_INSTANT(S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, Multiply_Epetra64, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_TPETRA_MATRIX_INSTANT)

#endif

}  // end namespace
