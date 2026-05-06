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
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Exceptions.hpp>

// #include <XpetraExt_MatrixMatrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractor.hpp>

#include <vector>

// #include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

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

RCP<const Comm<int>> getDefaultComm() {
  if (testMpi) {
    return DefaultPlatform::getDefaultPlatform().getComm();
  }
  return rcp(new Teuchos::SerialComm<int>());
}

/////////////////////////////////////////////////////

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
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

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, Multiply_Tpetra, M, MA, Scalar, LO, GO, Node) {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  // The matrix reader does not work with complex scalars. Skip all tests then.
  return;
#endif
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::CrsMatrix<Scalar, LO, GO, Node> CrsMatrixClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrapClass;

  // Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  // yAB->describe(*fos, Teuchos::VERB_EXTREME);

  {  // Tpetra test

    // get a comm and node
    RCP<const Comm<int>> comm = getDefaultComm();
    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    // define map
    LO nEle                       = 6;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    // read in matrices
    typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> reader_type;

    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpA    = reader_type::readSparseFile("A.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpB    = reader_type::readSparseFile("B.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpAB   = reader_type::readSparseFile("AB.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpAtB  = reader_type::readSparseFile("AtB.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpABt  = reader_type::readSparseFile("ABt.mat", comm);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpAtBt = reader_type::readSparseFile("AtBt.mat", comm);

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
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, BlockCrs, M, MB, Scalar, LO, GO, Node) {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  // The matrix reader does not work with complex scalars. Skip all tests then.
  return;
#endif
  typedef Tpetra::BlockCrsMatrix<Scalar, LO, GO, Node> BCM;
  typedef Tpetra::CrsGraph<LO, GO, Node> graph_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;

  typedef Xpetra::CrsMatrix<Scalar, LO, GO, Node> CrsMatrixClass;
  typedef Xpetra::TpetraBlockCrsMatrix<Scalar, LO, GO, Node> BlockCrsMatrixClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node> CrsMatrixWrapClass;
  using helpers = Xpetra::Helpers<Scalar, LO, GO, Node>;

  RCP<const Comm<int>> comm = getDefaultComm();
  const GO INVALID          = Teuchos::OrdinalTraits<GO>::invalid();

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
}

template <class Scalar, class LO, class GO, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>>
matrixAdd(const Xpetra::Matrix<Scalar, LO, GO, Node>& A,
          const Xpetra::Matrix<Scalar, LO, GO, Node>& B,
          Teuchos::FancyOStream& fos) {
  using MatrixClass           = Xpetra::Matrix<Scalar, LO, GO, Node>;
  Teuchos::RCP<MatrixClass> C = Teuchos::null;
  Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixAdd(
      A, false, Teuchos::ScalarTraits<Scalar>::one(),
      B, false, Teuchos::ScalarTraits<Scalar>::one(),
      C, fos);

  C->fillComplete(B.getDomainMap(), A.getRangeMap());
  return C;
}

template <class Scalar, class LO, class GO, class Node>
void comparePointMatrices(const Xpetra::Matrix<Scalar, LO, GO, Node>& expected,
                          const Xpetra::Matrix<Scalar, LO, GO, Node>& actual,
                          Teuchos::FancyOStream& out,
                          bool& success,
                          const std::string& label) {
  using ST      = Teuchos::ScalarTraits<Scalar>;
  using magType = typename ST::magnitudeType;
  using helpers = Xpetra::Helpers<Scalar, LO, GO, Node>;
  using crs_t   = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;

  TEST_ASSERT(helpers::isTpetraCrs(expected));
  TEST_ASSERT(helpers::isTpetraCrs(actual));

  const crs_t& exp = helpers::Op2TpetraCrs(expected);
  const crs_t& act = helpers::Op2TpetraCrs(actual);

  TEST_EQUALITY(exp.getLocalNumRows(), act.getLocalNumRows());
  TEST_EQUALITY(exp.getLocalNumEntries(), act.getLocalNumEntries());

  const magType tol = Teuchos::as<magType>(100) * Teuchos::ScalarTraits<magType>::eps();

  for (LO lrow = 0; lrow < static_cast<LO>(exp.getLocalNumRows()); ++lrow) {
    typename crs_t::local_inds_host_view_type expInds, actInds;
    typename crs_t::values_host_view_type expVals, actVals;

    exp.getLocalRowView(lrow, expInds, expVals);
    act.getLocalRowView(lrow, actInds, actVals);

    TEST_EQUALITY(expInds.extent(0), actInds.extent(0));

    for (size_t k = 0; k < expInds.extent(0); ++k) {
      TEST_EQUALITY(expInds[k], actInds[k]);
      TEST_FLOATING_EQUALITY(expVals[k], actVals[k], tol);
    }
  }

  out << "Verified " << label
      << " by local row comparison: local rows = " << exp.getLocalNumRows()
      << ", local nnz = " << exp.getLocalNumEntries()
      << std::endl;
}

template <class Scalar, class LO, class GO, class Node>
Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>>
build2x2BlockedMatrix(
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>>& A00,
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>>& A01,
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>>& A10,
    const Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>>& A11) {
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using map_type         = Xpetra::Map<LO, GO, Node>;
  using map_factory_type = Xpetra::MapFactory<LO, GO, Node>;
  using extractor_type   = Xpetra::MapExtractor<Scalar, LO, GO, Node>;
  using blocked_type     = Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>;

  TEUCHOS_TEST_FOR_EXCEPTION(A00.is_null(), Xpetra::Exceptions::RuntimeError,
                             "build2x2BlockedMatrix: A00 must be nonnull");

  RCP<const map_type> rowMap0 = A00->getRowMap();
  RCP<const map_type> domMap0 = A00->getDomainMap();

  auto lib       = rowMap0->lib();
  auto comm      = rowMap0->getComm();
  const GO rowIB = rowMap0->getIndexBase();
  const GO domIB = domMap0->getIndexBase();

  ArrayView<const GO> rowGids0 = rowMap0->getLocalElementList();
  ArrayView<const GO> domGids0 = domMap0->getLocalElementList();

  const GO rowOffset = rowMap0->getMaxAllGlobalIndex() + 1;
  const GO domOffset = domMap0->getMaxAllGlobalIndex() + 1;

  Array<GO> rowGids1(rowGids0.size()), domGids1(domGids0.size());
  for (size_t k = 0; k < static_cast<size_t>(rowGids0.size()); ++k)
    rowGids1[k] = rowGids0[k] + rowOffset;
  for (size_t k = 0; k < static_cast<size_t>(domGids0.size()); ++k)
    domGids1[k] = domGids0[k] + domOffset;

  RCP<const map_type> rowMap1 =
      map_factory_type::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                              rowGids1(), rowIB, comm);
  RCP<const map_type> domMap1 =
      map_factory_type::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                              domGids1(), domIB, comm);

  Array<GO> fullRowGids, fullDomGids;
  fullRowGids.reserve(rowGids0.size() + rowGids1.size());
  fullDomGids.reserve(domGids0.size() + domGids1.size());

  for (size_t k = 0; k < static_cast<size_t>(rowGids0.size()); ++k) fullRowGids.push_back(rowGids0[k]);
  for (size_t k = 0; k < static_cast<size_t>(rowGids1.size()); ++k) fullRowGids.push_back(rowGids1[k]);
  for (size_t k = 0; k < static_cast<size_t>(domGids0.size()); ++k) fullDomGids.push_back(domGids0[k]);
  for (size_t k = 0; k < static_cast<size_t>(domGids1.size()); ++k) fullDomGids.push_back(domGids1[k]);

  RCP<const map_type> fullRowMap =
      map_factory_type::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                              fullRowGids(), rowIB, comm);
  RCP<const map_type> fullDomMap =
      map_factory_type::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                              fullDomGids(), domIB, comm);

  std::vector<RCP<const map_type>> rowMaps(2), domMaps(2);
  rowMaps[0] = rowMap0;
  rowMaps[1] = rowMap1;
  domMaps[0] = domMap0;
  domMaps[1] = domMap1;

  RCP<const extractor_type> rangeExtractor =
      rcp(new extractor_type(fullRowMap, rowMaps, false));
  RCP<const extractor_type> domainExtractor =
      rcp(new extractor_type(fullDomMap, domMaps, false));

  RCP<blocked_type> blk = rcp(new blocked_type(rangeExtractor, domainExtractor, 33));

  blk->setMatrix(0, 0, A00);
  blk->setMatrix(0, 1, A01);
  blk->setMatrix(1, 0, A10);
  blk->setMatrix(1, 1, A11);
  blk->fillComplete();

  return blk;
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixMatrix, BlockedMultiply_Tpetra, M, MA, Scalar, LO, GO, Node) {
#if defined(HAVE_TPETRA_INST_COMPLEX_FLOAT) || defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  return;
#endif

  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using Xpetra::BlockedCrsMatrix;
  using Xpetra::CrsMatrix;
  using Xpetra::CrsMatrixWrap;
  using Xpetra::Matrix;

  using MatrixClass        = Matrix<Scalar, LO, GO, Node>;
  using CrsMatrixClass     = CrsMatrix<Scalar, LO, GO, Node>;
  using CrsMatrixWrapClass = CrsMatrixWrap<Scalar, LO, GO, Node>;
  using BlockedClass       = BlockedCrsMatrix<Scalar, LO, GO, Node>;

  RCP<const Comm<int>> comm = getDefaultComm();

  typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> reader_type;

  RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpA    = reader_type::readSparseFile("A.mat", comm);
  RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpB    = reader_type::readSparseFile("B.mat", comm);
  RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpAB   = reader_type::readSparseFile("AB.mat", comm);
  RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpAtB  = reader_type::readSparseFile("AtB.mat", comm);
  RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpABt  = reader_type::readSparseFile("ABt.mat", comm);
  RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node>> tpAtBt = reader_type::readSparseFile("AtBt.mat", comm);

  RCP<CrsMatrixClass> xAmat    = rcp(new MA(tpA));
  RCP<CrsMatrixClass> xBmat    = rcp(new MA(tpB));
  RCP<CrsMatrixClass> xABmat   = rcp(new MA(tpAB));
  RCP<CrsMatrixClass> xAtBmat  = rcp(new MA(tpAtB));
  RCP<CrsMatrixClass> xABtmat  = rcp(new MA(tpABt));
  RCP<CrsMatrixClass> xAtBtmat = rcp(new MA(tpAtBt));

  RCP<MatrixClass> xA    = rcp(new CrsMatrixWrapClass(xAmat));
  RCP<MatrixClass> xB    = rcp(new CrsMatrixWrapClass(xBmat));
  RCP<MatrixClass> xAB   = rcp(new CrsMatrixWrapClass(xABmat));
  RCP<MatrixClass> xAtB  = rcp(new CrsMatrixWrapClass(xAtBmat));
  RCP<MatrixClass> xABt  = rcp(new CrsMatrixWrapClass(xABtmat));
  RCP<MatrixClass> xAtBt = rcp(new CrsMatrixWrapClass(xAtBtmat));

  auto fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(out));

  //
  // Ablk = [ A  A ]
  //        [ 0  A ]
  //
  // Bblk = [ B  0 ]
  //        [ B  B ]
  //
  RCP<BlockedClass> Ablk = build2x2BlockedMatrix<Scalar, LO, GO, Node>(xA, xA, Teuchos::null, xA);
  RCP<BlockedClass> Bblk = build2x2BlockedMatrix<Scalar, LO, GO, Node>(xB, Teuchos::null, xB, xB);

  RCP<MatrixClass> twoAB   = matrixAdd<Scalar, LO, GO, Node>(*xAB, *xAB, *fos);
  RCP<MatrixClass> twoAtB  = matrixAdd<Scalar, LO, GO, Node>(*xAtB, *xAtB, *fos);
  RCP<MatrixClass> twoABt  = matrixAdd<Scalar, LO, GO, Node>(*xABt, *xABt, *fos);
  RCP<MatrixClass> twoAtBt = matrixAdd<Scalar, LO, GO, Node>(*xAtBt, *xAtBt, *fos);

  struct CaseData {
    bool transposeA, transposeB;
    RCP<MatrixClass> e00, e01, e10, e11;
    bool has00, has01, has10, has11;
    std::string label;
  };

  Teuchos::Array<CaseData> cases;
  cases.push_back({false, false, twoAB, xAB, xAB, xAB, true, true, true, true, "NN"});
  cases.push_back({true, false, xAtB, Teuchos::null, twoAtB, xAtB, true, false, true, true, "TN"});
  cases.push_back({false, true, xABt, twoABt, Teuchos::null, xABt, true, true, false, true, "NT"});
  cases.push_back({true, true, xAtBt, xAtBt, xAtBt, twoAtBt, true, true, true, true, "TT"});

  for (const auto& tc : cases) {
    out << "Testing blocked multiply case " << tc.label
        << " transposeA=" << tc.transposeA
        << " transposeB=" << tc.transposeB << std::endl;

    RCP<BlockedClass> Cblk =
        Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixMultiplyBlock(
            *Ablk, tc.transposeA, *Bblk, tc.transposeB, *fos, true, true);

    TEST_ASSERT(!Cblk.is_null());
    TEST_EQUALITY(Cblk->Rows(), static_cast<size_t>(2));
    TEST_EQUALITY(Cblk->Cols(), static_cast<size_t>(2));

    RCP<MatrixClass> C00 = Cblk->getMatrix(0, 0);
    RCP<MatrixClass> C01 = Cblk->getMatrix(0, 1);
    RCP<MatrixClass> C10 = Cblk->getMatrix(1, 0);
    RCP<MatrixClass> C11 = Cblk->getMatrix(1, 1);

    TEST_EQUALITY(!C00.is_null(), tc.has00);
    TEST_EQUALITY(!C01.is_null(), tc.has01);
    TEST_EQUALITY(!C10.is_null(), tc.has10);
    TEST_EQUALITY(!C11.is_null(), tc.has11);

    if (tc.has00) comparePointMatrices<Scalar, LO, GO, Node>(*tc.e00, *C00, out, success, tc.label + "(0,0)");
    if (tc.has01) comparePointMatrices<Scalar, LO, GO, Node>(*tc.e01, *C01, out, success, tc.label + "(0,1)");
    if (tc.has10) comparePointMatrices<Scalar, LO, GO, Node>(*tc.e10, *C10, out, success, tc.label + "(1,0)");
    if (tc.has11) comparePointMatrices<Scalar, LO, GO, Node>(*tc.e11, *C11, out, success, tc.label + "(1,1)");
  }
}

//
// INSTANTIATIONS
//

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                                  \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N;              \
  typedef typename Xpetra::TpetraCrsMatrix<S, LO, GO, N> MA##S##LO##GO##N; \
  typedef typename Xpetra::TpetraBlockCrsMatrix<S, LO, GO, N> MB##S##LO##GO##N;

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S, LO, GO, N)                                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, Multiply_Tpetra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, BlockCrs, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixMatrix, BlockedMultiply_Tpetra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_TPETRA_MATRIX_INSTANT)

}  // end namespace
