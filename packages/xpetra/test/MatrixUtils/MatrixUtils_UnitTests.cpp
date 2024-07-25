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

#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MatrixUtils.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsMatrix.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif
namespace {
using Xpetra::viewLabel_t;

bool testMpi         = true;
double errorTolSlack = 1e+1;

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

Teuchos::RCP<const Teuchos::Comm<int>> getDefaultComm() {
  if (testMpi) {
    return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  }
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
}

//
// UNIT TESTS
//

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MatrixUtils, SwitchMatrixToStridedMaps, M, MA, Scalar, LO, GO, Node) {
  using CrsMatrixWrap     = Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>;
  using Map               = Xpetra::Map<LO, GO, Node>;
  using MapFactory        = Xpetra::MapFactory<LO, GO, Node>;
  using MatrixUtils       = Xpetra::MatrixUtils<Scalar, LO, GO, Node>;
  using StridedMap        = Xpetra::StridedMap<LO, GO, Node>;
  using StridedMapFactory = Xpetra::StridedMapFactory<LO, GO, Node>;

  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();
  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const size_t numLocalElements = 10;
  const GO INVALID              = Teuchos::OrdinalTraits<GO>::invalid();
  RCP<const Map> map            = MapFactory::createContigMapWithNode(lib, INVALID, numLocalElements, comm);

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(1);
  RCP<const StridedMap> stridedMap = StridedMapFactory::Build(map, stridingInfo);

  RCP<CrsMatrixWrap> matrix = rcp(new CrsMatrixWrap(stridedMap, stridedMap, 1));
  TEST_EQUALITY_CONST(matrix->GetCurrentViewLabel(), matrix->GetDefaultViewLabel());
  TEST_EQUALITY_CONST(matrix->GetCurrentViewLabel(), matrix->SwitchToView(matrix->GetCurrentViewLabel()));

  for (LO lRow = 0; lRow < Teuchos::as<LO>(numLocalElements); ++lRow) {
    Array<LO> lCols(1, lRow);
    Array<Scalar> values(1, Teuchos::ScalarTraits<Scalar>::one());
    matrix->insertLocalValues(lRow, lCols, values);
  }

  matrix->fillComplete();
  TEST_ASSERT(matrix->isFillComplete());

  TEST_EQUALITY_CONST(matrix->IsView("stridedMaps"), false)
  MatrixUtils::convertMatrixToStridedMaps(matrix, stridingInfo, stridingInfo);
  TEST_EQUALITY_CONST(matrix->IsView("stridedMaps"), true)

  TEST_ASSERT(!matrix->getRangeMap().is_null());
  TEST_ASSERT(!matrix->getDomainMap().is_null());
  TEST_ASSERT(!matrix->getRowMap().is_null());
  TEST_ASSERT(!matrix->getColMap().is_null());

  matrix->SwitchToView("stridedMaps");
  RCP<const StridedMap> stridedRowMap = rcp_dynamic_cast<const StridedMap>(matrix->getRowMap("stridedMaps"));
  RCP<const StridedMap> stridedColMap = rcp_dynamic_cast<const StridedMap>(matrix->getColMap("stridedMaps"));

  TEST_ASSERT(!stridedRowMap.is_null());
  TEST_ASSERT(!stridedColMap.is_null());
}

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                     \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N; \
  typedef typename Xpetra::TpetraCrsMatrix<S, LO, GO, N> MA##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

#define XPETRA_EPETRA_TYPES(S, LO, GO, N)                  \
  typedef typename Xpetra::EpetraMapT<GO, N> M##LO##GO##N; \
  typedef typename Xpetra::EpetraCrsMatrixT<GO, N> MA##S##LO##GO##N;

#endif

// list of all tests which run both with Epetra and Tpetra
#define XP_MATRIX_INSTANT(S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MatrixUtils, SwitchMatrixToStridedMaps, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_MATRIX_INSTANT)

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double, int, int, EpetraNode)
XP_MATRIX_INSTANT(double, int, int, EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double, int, LongLong, EpetraNode)
XP_MATRIX_INSTANT(double, int, LongLong, EpetraNode)
#endif

#endif

}  // namespace
