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
#include <Teuchos_as.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Tuple.hpp>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Epetra_SerialComm.h"

#include <Xpetra_ConfigDefs.hpp>

#ifdef HAVE_XPETRA_EPETRAEXT
// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

// Epetra routines to split matrix and maps
#include "BlockedMatrixTestHelpers.hpp"

#include <Xpetra_DefaultPlatform.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapUtils.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_ReorderedBlockedMultiVector.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

namespace XpetraBlockMatrixTests {

bool testMpi         = true;
double errorTolSlack = 1e+1;

Teuchos::RCP<const Teuchos::Comm<int>> getDefaultComm() {
  if (testMpi) {
    return Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
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

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, SplitMatrix, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar, LO, GO, Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar, LO, GO, Node> MapExtractorFactoryClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  GO nEle                                = 63;
  const Teuchos::RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  LO NumMyElements                              = map->getLocalNumElements();
  GO NumGlobalElements                          = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (LO i = 0; i < NumMyElements; i++) {
    if (MyGlobalElements[i] == 0) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(Teuchos::as<Scalar>(i) * STS::one(), -1.0));
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(-1.0, Teuchos::as<Scalar>(i) * STS::one()));
    } else {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(-1.0, Teuchos::as<Scalar>(i) * STS::one(), -1.0));
    }
  }

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> mat =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>(A));

  Teuchos::Array<GO> gids1;
  Teuchos::Array<GO> gids2;
  for (LO i = 0; i < NumMyElements; i++) {
    if (i % 3 < 2)
      gids1.push_back(map->getGlobalElement(i));
    else
      gids2.push_back(map->getGlobalElement(i));
  }

  const Teuchos::RCP<const MapClass> map1 = MapFactoryClass::Build(lib,
                                                                   Teuchos::OrdinalTraits<GO>::invalid(),
                                                                   gids1.view(0, gids1.size()),
                                                                   0,
                                                                   comm);
  const Teuchos::RCP<const MapClass> map2 = MapFactoryClass::Build(lib,
                                                                   Teuchos::OrdinalTraits<GO>::invalid(),
                                                                   gids2.view(0, gids2.size()),
                                                                   0,
                                                                   comm);

  std::vector<Teuchos::RCP<const MapClass>> xmaps;
  xmaps.push_back(map1);
  xmaps.push_back(map2);

  Teuchos::RCP<const MapExtractorClass> rowMapExtractormap_extractor = MapExtractorFactoryClass::Build(map, xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bOp =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*mat, rowMapExtractormap_extractor, rowMapExtractormap_extractor);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(map, true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(map, true);
  ones->putScalar(STS::one());
  rnd->randomize();

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  magnitudeType tol1  = Teuchos::ScalarTraits<magnitudeType>::eps();
  magnitudeType tol2  = 500 * tol1;

  A->apply(*ones, *exp);
  bOp->apply(*ones, *res);
  res->update(-STS::one(), *exp, STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, tol1, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, tol1, out, success);

  A->apply(*rnd, *exp);
  bOp->apply(*rnd, *res);
  res->update(-STS::one(), *exp, STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, tol2, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, tol2, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, CreateBlockedDiagonalOp, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 4;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, noBlocks - 2)) * 10 * comm->getSize();

  TEST_EQUALITY(bop->Rows(), 4);
  TEST_EQUALITY(bop->Cols(), 4);
  TEST_EQUALITY(bop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(bop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 9);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40 + 20);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getLocalNumElements(), 10);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getLocalNumElements(), 20);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getMinGlobalIndex(), comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getMinGlobalIndex(), comm->getRank() * 40 + 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 9);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getMinGlobalIndex(), comm->getRank() * 40 + 20);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getLocalNumElements(), 10);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getLocalNumElements(), 20);

  TEST_EQUALITY(bop->getMatrix(0, 1)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0, 1)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(0, 2)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0, 2)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(0, 3)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(0, 3)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  TEST_EQUALITY(bop->getMatrix(1, 0)->getColMap()->getMinGlobalIndex(), std::numeric_limits<GO>::max());  // TODO
  TEST_EQUALITY(bop->getMatrix(1, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  // TEST_EQUALITY(bop->getMatrix(1,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(2, 0)->getColMap()->getMinGlobalIndex(), std::numeric_limits<GO>::max());  // TODO
  TEST_EQUALITY(bop->getMatrix(2, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  // TEST_EQUALITY(bop->getMatrix(2,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40);
  TEST_EQUALITY(bop->getMatrix(3, 0)->getColMap()->getMinGlobalIndex(), std::numeric_limits<GO>::max());  // TODO
  TEST_EQUALITY(bop->getMatrix(3, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 4);
  // TEST_EQUALITY(bop->getMatrix(3,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 40);

  TEST_EQUALITY(bop->getMatrix(2, 1)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2, 1)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);
  TEST_EQUALITY(bop->getMatrix(2, 3)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 40 + 10);
  TEST_EQUALITY(bop->getMatrix(2, 3)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 19);

  TEST_EQUALITY(bop->getMatrix(0, 0)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(2, 2)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(2, 3)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(1, 0)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(3, 1)->isFillComplete(), true);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bop->getDomainMapExtractor()->getThyraMode(), false);

#ifdef HAVE_XPETRA_DEBUG
  TEST_THROW(bop->getRangeMap(0, true), Xpetra::Exceptions::RuntimeError);
#endif

  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(2)->getMinAllGlobalIndex(), 10);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(2)->getMaxAllGlobalIndex(), (comm->getSize() - 1) * 40 + 19);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMinGlobalIndex(), comm->getRank() * 40 + 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMinAllGlobalIndex(), 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, CreateBlockedDiagonalOpThyra, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 4;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, noBlocks - 2)) * 10 * comm->getSize();

  TEST_EQUALITY(bop->Rows(), 4);
  TEST_EQUALITY(bop->Cols(), 4);
  TEST_EQUALITY(bop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(bop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  // Thyra GIDs
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 20);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 20 + 19);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(bop->getMatrix(0, 0)->getRowMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getRowMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getLocalNumElements(), 10);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getRowMap()->getLocalNumElements(), 20);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getMinGlobalIndex(), comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getMinGlobalIndex(), comm->getRank() * 20);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getMaxGlobalIndex(), comm->getRank() * 20 + 19);
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(bop->getMatrix(0, 0)->getColMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(1, 1)->getColMap()->getLocalNumElements(), 5);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getColMap()->getLocalNumElements(), 10);
  TEST_EQUALITY(bop->getMatrix(3, 3)->getColMap()->getLocalNumElements(), 20);

  TEST_EQUALITY(bop->getMatrix(0, 1)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0, 1)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(0, 2)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0, 2)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(0, 3)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(0, 3)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getMatrix(1, 0)->getColMap()->getMinGlobalIndex(), std::numeric_limits<GO>::max());  // TODO
  TEST_EQUALITY(bop->getMatrix(1, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  // TEST_EQUALITY(bop->getMatrix(1,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(2, 0)->getColMap()->getMinGlobalIndex(), std::numeric_limits<GO>::max());  // TODO
  TEST_EQUALITY(bop->getMatrix(2, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  // TEST_EQUALITY(bop->getMatrix(2,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5);
  TEST_EQUALITY(bop->getMatrix(3, 0)->getColMap()->getMinGlobalIndex(), std::numeric_limits<GO>::max());  // TODO
  TEST_EQUALITY(bop->getMatrix(3, 0)->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  // TEST_EQUALITY(bop->getMatrix(3,0)->getColMap()->getMaxGlobalIndex(),comm->getRank() * 5);

  TEST_EQUALITY(bop->getMatrix(2, 1)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2, 1)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2, 2)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getMatrix(2, 3)->getRowMap()->getMinGlobalIndex(), comm->getRank() * 10);
  TEST_EQUALITY(bop->getMatrix(2, 3)->getRowMap()->getMaxGlobalIndex(), comm->getRank() * 10 + 9);

  TEST_EQUALITY(bop->getMatrix(0, 0)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(2, 2)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(2, 3)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(1, 0)->isFillComplete(), true);
  TEST_EQUALITY(bop->getMatrix(3, 1)->isFillComplete(), true);

  // check Xpetra replacement maps
  TEST_EQUALITY(bop->getRangeMap(0, false)->getMinGlobalIndex(), comm->getRank() * 5 + 0);
  TEST_EQUALITY(bop->getRangeMap(0, false)->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(1, false)->getMinGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5);
  TEST_EQUALITY(bop->getRangeMap(1, false)->getMaxGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(2, false)->getMinGlobalIndex(), comm->getSize() * 10 + comm->getRank() * 10);
  TEST_EQUALITY(bop->getRangeMap(2, false)->getMaxGlobalIndex(), comm->getSize() * 10 + comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getRangeMap(3, false)->getMinGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(bop->getRangeMap(3, false)->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);

  // check Thyra maps
  TEST_EQUALITY(bop->getRangeMap(0)->getMinGlobalIndex(), comm->getRank() * 5 + 0);
  TEST_EQUALITY(bop->getRangeMap(0)->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(1)->getMinGlobalIndex(), comm->getRank() * 5 + 0);
  TEST_EQUALITY(bop->getRangeMap(1)->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(bop->getRangeMap(2)->getMinGlobalIndex(), comm->getRank() * 10 + 0);
  TEST_EQUALITY(bop->getRangeMap(2)->getMaxGlobalIndex(), comm->getRank() * 10 + 9);
  TEST_EQUALITY(bop->getRangeMap(3)->getMinGlobalIndex(), comm->getRank() * 20 + 0);
  TEST_EQUALITY(bop->getRangeMap(3)->getMaxGlobalIndex(), comm->getRank() * 20 + 19);

  TEST_EQUALITY(bop->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bop->getRangeMap(0)->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);
  TEST_EQUALITY(bop->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bop->getRangeMap(1)->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);
  TEST_EQUALITY(bop->getRangeMap(2)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bop->getRangeMap(2)->getMaxAllGlobalIndex(), comm->getSize() * 10 - 1);
  TEST_EQUALITY(bop->getRangeMap(3)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bop->getRangeMap(3)->getMaxAllGlobalIndex(), comm->getSize() * 20 - 1);

  TEST_EQUALITY(bop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bop->getDomainMapExtractor()->getThyraMode(), true);

  // check Xpetra replacement submaps
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMinGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMinAllGlobalIndex(), comm->getSize() * 20);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getMap(3, false)->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperator, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 0 [ [1 2] 3] ] 4 [ 5 6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, noBlocks - 2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));

  // block 00
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop00 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop->getMatrix(0, 0));

  GO goNumRows00 = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, 2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop00->Rows(), 2);
  TEST_EQUALITY(brop00->Cols(), 2);
  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows00));
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows00));

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop11 = brop->getMatrix(1, 1);

  GO goNumRows11 = Teuchos::as<GO>(40 * comm->getSize());
  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows11));
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows11));
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 79);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop11test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop11);
  TEST_INEQUALITY(brop11test, Teuchos::null);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  TEST_EQUALITY(brop11test->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows11));
  TEST_EQUALITY(brop11test->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows11));
  TEST_EQUALITY(brop11test->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop11test->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 79);

  // block 22
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop22 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop->getMatrix(2, 2));

  GO goNumRows22 = Teuchos::as<GO>(560 * comm->getSize());

  TEST_EQUALITY(brop22->Rows(), 3);
  TEST_EQUALITY(brop22->Cols(), 3);
  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows22));
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows22));
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 80);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 639);
  TEST_EQUALITY(brop22->getMatrix(0, 0)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 80);
  TEST_EQUALITY(brop22->getMatrix(0, 0)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 159);
  TEST_EQUALITY(brop22->getMatrix(1, 1)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 160);
  TEST_EQUALITY(brop22->getMatrix(1, 1)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 319);
  TEST_EQUALITY(brop22->getMatrix(2, 2)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 320);
  TEST_EQUALITY(brop22->getMatrix(2, 2)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 639);

  // block 00_11
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop00_11 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop00->getMatrix(1, 1));

  GO goNumRows00_11 = Teuchos::as<GO>(35 * comm->getSize());

  TEST_EQUALITY(brop00_11->Rows(), 2);
  TEST_EQUALITY(brop00_11->Cols(), 2);
  TEST_EQUALITY(brop00_11->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows00_11));
  TEST_EQUALITY(brop00_11->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows00_11));
  TEST_EQUALITY(brop00_11->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 39);
  TEST_EQUALITY(brop00_11->getMatrix(0, 0)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(15 * comm->getSize()));
  TEST_EQUALITY(brop00_11->getMatrix(0, 0)->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(15 * comm->getSize()));
  TEST_EQUALITY(brop00_11->getMatrix(0, 0)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getMatrix(0, 0)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 19);
  TEST_EQUALITY(brop00_11->getMatrix(1, 1)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(20 * comm->getSize()));
  TEST_EQUALITY(brop00_11->getMatrix(1, 1)->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(20 * comm->getSize()));
  TEST_EQUALITY(brop00_11->getMatrix(1, 1)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop00_11->getMatrix(1, 1)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 39);

  // block 01
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop01 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop->getMatrix(0, 1));

  TEST_EQUALITY(brop01->Rows(), 2);
  TEST_EQUALITY(brop01->Cols(), 1);
  TEST_EQUALITY(brop01->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop01->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop01->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 0);
  TEST_EQUALITY(brop01->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 39);
  TEST_EQUALITY(brop01->getDomainMap()->getMinGlobalIndex(), comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop01->getDomainMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 79);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperator2, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 4 3 1 7 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 4);
  TEST_EQUALITY(brop->Cols(), 4);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 385));
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 385));

  // block 00
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop00 = brop->getMatrix(0, 0);

  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop00->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 79);

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop11 = brop->getMatrix(1, 1);

  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 39);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop11test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop11);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  TEST_EQUALITY(brop11test->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(brop11test->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(brop11test->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop11test->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 39);

  // block 22
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop22 = brop->getMatrix(2, 2);

  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 9);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop22test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop22);

  TEST_EQUALITY(brop22test->Rows(), 1);
  TEST_EQUALITY(brop22test->Cols(), 1);
  TEST_EQUALITY(brop22test->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(brop22test->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(brop22test->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop22test->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 9);

  // block 33
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop33 = brop->getMatrix(3, 3);

  TEST_EQUALITY(brop33->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 320));
  TEST_EQUALITY(brop33->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 320));
  TEST_EQUALITY(brop33->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 320);
  TEST_EQUALITY(brop33->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 639);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop33test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop33);

  TEST_EQUALITY(brop33test->Rows(), 1);
  TEST_EQUALITY(brop33test->Cols(), 1);
  TEST_EQUALITY(brop33test->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 320));
  TEST_EQUALITY(brop33test->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 320));
  TEST_EQUALITY(brop33test->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 640 + 320);
  TEST_EQUALITY(brop33test->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 640 + 639);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperatorThyra, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> Matrix;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node> ReorderedBlockedCrsMatrix;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 0 [ [1 2] 3] ] 4 [ 5 6 7] ]");

  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(buildReorderedBlockedCrsMatrix(brm, bop));

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, noBlocks - 2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));

  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // block 00
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop00 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop->getMatrix(0, 0));

  GO goNumRows00 = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, 2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop00->Rows(), 2);
  TEST_EQUALITY(brop00->Cols(), 2);
  TEST_EQUALITY(brop00->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop00->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows00));
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows00));
  TEST_EQUALITY(brop00->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 5 + 0);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop00->getRangeMap()->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);
  // Thyra maps (since it is a blocked matrix, they should be unique!)
  TEST_EQUALITY(brop00->getRangeMap(0, true)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap(0, true)->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);
  TEST_EQUALITY(brop00->getRangeMap(1, true)->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 35));
  TEST_EQUALITY(brop00->getRangeMap(1, true)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap(1, true)->getMaxAllGlobalIndex(), comm->getSize() * 35 - 1);
  // Xpetra maps
  TEST_EQUALITY(brop00->getRangeMap(0, false)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop00->getRangeMap(0, false)->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);
  TEST_EQUALITY(brop00->getRangeMap(1, false)->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 35));
  TEST_EQUALITY(brop00->getRangeMap(1, false)->getMinAllGlobalIndex(), comm->getSize() * 5);
  TEST_EQUALITY(brop00->getRangeMap(1, false)->getMaxAllGlobalIndex(), comm->getSize() * 5 + comm->getSize() * 35 - 1);

  // subblock 11 of block 00
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> sbrop11 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop00->getMatrix(1, 1));

  TEST_EQUALITY(sbrop11->Rows(), 2);
  TEST_EQUALITY(sbrop11->Cols(), 2);
  TEST_EQUALITY(sbrop11->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(sbrop11->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(sbrop11->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(35 * comm->getSize()));
  TEST_EQUALITY(sbrop11->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(35 * comm->getSize()));
  TEST_EQUALITY(sbrop11->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);
  // Thyra maps (since it is a blocked matrix, they should be unique!)
  TEST_EQUALITY(sbrop11->getRangeMap(0, true)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(sbrop11->getRangeMap(0, true)->getMaxAllGlobalIndex(), comm->getSize() * 5 + comm->getSize() * 10 - 1);
  TEST_EQUALITY(sbrop11->getRangeMap(1, true)->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(sbrop11->getRangeMap(1, true)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(sbrop11->getRangeMap(1, true)->getMaxAllGlobalIndex(), comm->getSize() * 20 - 1);
  // Xpetra maps
  TEST_EQUALITY(sbrop11->getRangeMap(0, false)->getMinAllGlobalIndex(), comm->getSize() * 5);
  TEST_EQUALITY(sbrop11->getRangeMap(0, false)->getMaxAllGlobalIndex(), comm->getSize() * 20 - 1);
  TEST_EQUALITY(sbrop11->getRangeMap(1, false)->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(sbrop11->getRangeMap(1, false)->getMinAllGlobalIndex(), comm->getSize() * 20);
  TEST_EQUALITY(sbrop11->getRangeMap(1, false)->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop11 = brop->getMatrix(1, 1);

  // Thyra GIDs for the matrix
  GO goNumRows11 = Teuchos::as<GO>(40 * comm->getSize());
  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows11));
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows11));
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 0);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop11->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 40 + 0);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 + comm->getSize() * 40 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getRangeMap(1, false)->getMinAllGlobalIndex(), comm->getSize() * 40);
  TEST_EQUALITY(brop->getRangeMap(1, false)->getMaxAllGlobalIndex(), 2 * comm->getSize() * 40 - 1);

  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop11test =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop11);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  TEST_EQUALITY(brop11test->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop11test->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 40 + 0);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 - 1);

  Teuchos::RCP<Matrix> crsmat11 = Teuchos::rcp_const_cast<ReorderedBlockedCrsMatrix>(brop11test)->getInnermostCrsMatrix();
  TEST_EQUALITY(crsmat11.is_null(), false);
  Teuchos::ArrayView<const LO> inds;
  Teuchos::ArrayView<const Scalar> vals;
  crsmat11->getLocalRowView(0, inds, vals);
  TEST_EQUALITY(inds.size(), 1);
  TEST_EQUALITY(vals.size(), 1);
  TEST_EQUALITY(inds[0], 0);
  TEST_EQUALITY(vals[0], Teuchos::as<Scalar>(5.0) * Teuchos::ScalarTraits<Scalar>::one());

  // block 22
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop22 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop->getMatrix(2, 2));

  GO goNumRows22 = Teuchos::as<GO>(560 * comm->getSize());

  TEST_EQUALITY(brop22->Rows(), 3);
  TEST_EQUALITY(brop22->Cols(), 3);
  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows22));
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows22));
  // Xpetra replacement GIDs
  TEST_EQUALITY(brop22->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 80);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 80 + comm->getSize() * 560 - 1);
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 80 + comm->getSize() * 80);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 80 + comm->getSize() * 240 + comm->getRank() * 320 + 319);
  // Xpetra GIDs
  TEST_EQUALITY(brop22->getRangeMap(0, false)->getMinGlobalIndex(), comm->getSize() * 80 + comm->getRank() * 80);
  TEST_EQUALITY(brop22->getRangeMap(0, false)->getMaxGlobalIndex(), comm->getSize() * 80 + comm->getRank() * 80 + 79);
  TEST_EQUALITY(brop22->getRangeMap(1, false)->getMinGlobalIndex(), comm->getSize() * 160 + comm->getRank() * 160);
  TEST_EQUALITY(brop22->getRangeMap(1, false)->getMaxGlobalIndex(), comm->getSize() * 160 + comm->getRank() * 160 + 159);
  TEST_EQUALITY(brop22->getRangeMap(2, false)->getMinGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320);
  TEST_EQUALITY(brop22->getRangeMap(2, false)->getMaxGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320 + 319);

  // block 00_11
  /*Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> > brop00_11 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar,LO,GO,Node> >(brop00->getMatrix(1,1));

  GO goNumRows00_11 = Teuchos::as<GO>(35 * comm->getSize());

  TEST_EQUALITY(brop00_11->Rows(),2);
  TEST_EQUALITY(brop00_11->Cols(),2);
  TEST_EQUALITY(brop00_11->getRangeMap()->getGlobalNumElements(),goNumRows00_11);
  TEST_EQUALITY(brop00_11->getDomainMap()->getGlobalNumElements(),goNumRows00_11);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getGlobalNumElements(),15 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getDomainMap()->getGlobalNumElements(),15 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 5);
  TEST_EQUALITY(brop00_11->getMatrix(0,0)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 19);
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getGlobalNumElements(),20 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getDomainMap()->getGlobalNumElements(),20 * comm->getSize());
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getMinGlobalIndex(),comm->getRank() * 640 + 20);
  TEST_EQUALITY(brop00_11->getMatrix(1,1)->getRangeMap()->getMaxGlobalIndex(),comm->getRank() * 640 + 39);
*/
  // block 01
  Teuchos::RCP<const ReorderedBlockedCrsMatrix> brop01 =
      Teuchos::rcp_dynamic_cast<const ReorderedBlockedCrsMatrix>(brop->getMatrix(0, 1));

  // Xpetra like maps
  TEST_EQUALITY(brop01->Rows(), 2);
  TEST_EQUALITY(brop01->Cols(), 1);
  TEST_EQUALITY(brop01->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop01->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop01->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 5 + 0);
  TEST_EQUALITY(brop01->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 5 + comm->getSize() * 15 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop01->getDomainMap()->getMinGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop01->getDomainMap()->getMaxGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 39);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperator2Thyra, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 4 3 1 7 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 4);
  TEST_EQUALITY(brop->Cols(), 4);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 385));
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 385));

  // block 00
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop00 = brop->getMatrix(0, 0);

  TEST_EQUALITY(brop00->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop00->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  // Thyra GIDs
  TEST_EQUALITY(brop00->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 39);
  TEST_EQUALITY(brop00->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 40);
  TEST_EQUALITY(brop00->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 40 + comm->getSize() * 40 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(0, false)->getMinGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40);
  TEST_EQUALITY(brop->getDomainMap(0, false)->getMaxGlobalIndex(), comm->getSize() * 40 + comm->getRank() * 40 + 39);

  // block 11
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop11 = brop->getMatrix(1, 1);

  TEST_EQUALITY(brop11->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  TEST_EQUALITY(brop11->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 20));
  // Thyra GIDs (+ Xpetra shift)
  TEST_EQUALITY(brop11->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop11->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 20 + 0);
  TEST_EQUALITY(brop11->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 20 + comm->getSize() * 20 - 1);

  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(1, false)->getMinGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20);
  TEST_EQUALITY(brop->getDomainMap(1, false)->getMaxGlobalIndex(), comm->getSize() * 20 + comm->getRank() * 20 + 19);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop11test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop11);

  TEST_EQUALITY(brop11test->Rows(), 1);
  TEST_EQUALITY(brop11test->Cols(), 1);
  // Thyra GIDs
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 20);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 20 + 19);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop11test->getMatrix(0, 0)->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 20 - 1);

  // block 22
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop22 = brop->getMatrix(2, 2);

  TEST_EQUALITY(brop22->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(brop22->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));
  TEST_EQUALITY(brop22->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5 + 4);
  TEST_EQUALITY(brop22->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 5 + 0);
  TEST_EQUALITY(brop22->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 5 + comm->getSize() * 5 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(2, false)->getMinGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5);
  TEST_EQUALITY(brop->getDomainMap(2, false)->getMaxGlobalIndex(), comm->getSize() * 5 + comm->getRank() * 5 + 4);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop22test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop22);
  TEST_EQUALITY(brop22test->Rows(), 1);
  TEST_EQUALITY(brop22test->Cols(), 1);

  // Thyra GIDs
  TEST_EQUALITY(brop22test->getMatrix(0, 0)->getRangeMap()->getMinGlobalIndex(), comm->getRank() * 5);
  TEST_EQUALITY(brop22test->getMatrix(0, 0)->getRangeMap()->getMaxGlobalIndex(), comm->getRank() * 5 + 4);
  TEST_EQUALITY(brop22test->getMatrix(0, 0)->getRangeMap()->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(brop22test->getMatrix(0, 0)->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 5 - 1);

  // block 33
  Teuchos::RCP<const Xpetra::Matrix<Scalar, LO, GO, Node>> brop33 = brop->getMatrix(3, 3);

  TEST_EQUALITY(brop33->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 320));
  TEST_EQUALITY(brop33->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 320));
  TEST_EQUALITY(brop33->getRangeMap()->getMinGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320);
  TEST_EQUALITY(brop33->getRangeMap()->getMaxGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320 + 319);
  TEST_EQUALITY(brop33->getRangeMap()->getMinAllGlobalIndex(), comm->getSize() * 320 + 0);
  TEST_EQUALITY(brop33->getRangeMap()->getMaxAllGlobalIndex(), comm->getSize() * 320 + comm->getSize() * 320 - 1);
  // Xpetra GIDs
  TEST_EQUALITY(brop->getDomainMap(3, false)->getMinGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320);
  TEST_EQUALITY(brop->getDomainMap(3, false)->getMaxGlobalIndex(), comm->getSize() * 320 + comm->getRank() * 320 + 319);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop33test =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop33);
  TEST_EQUALITY(brop33test->Rows(), 1);
  TEST_EQUALITY(brop33test->Cols(), 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperatorApply, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVectorClass;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 [ 2 3 4 ] 5 ] [6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);

  // build gloabl vector with one entries
  // The MultiVector objects "ones" and "exp" are BlockedMultiVectors with 8 sub blocks
  // compatible to bop
  Teuchos::RCP<MultiVectorClass> ones = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  Teuchos::RCP<MultiVectorClass> exp  = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  ones->putScalar(STS::one());
  bop->apply(*ones, *exp);

  // reorganize "ones" and "res" to be BlockedMultiVectors with 3 sub blocks (nested)
  // compatible to brop
  // They use and work with the same 8 sub vectors from "ones" and "exp"
  Teuchos::RCP<const MultiVectorClass> cones = Teuchos::rcp_const_cast<const MultiVectorClass>(ones);
  Teuchos::RCP<const MultiVectorClass> brones =
      Xpetra::buildReorderedBlockedMultiVector(brm, Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(cones));
  Teuchos::RCP<MultiVectorClass> res        = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  Teuchos::RCP<const MultiVectorClass> cres = Teuchos::rcp_const_cast<const MultiVectorClass>(res);
  Teuchos::RCP<const MultiVectorClass> brcres =
      Xpetra::buildReorderedBlockedMultiVector(brm, Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(cres));
  Teuchos::RCP<MultiVectorClass> brres = Teuchos::rcp_const_cast<MultiVectorClass>(brcres);

  brop->apply(*brones, *brres);

  Teuchos::Array<typename STS::magnitudeType> nn(res->getNumVectors());
  TEST_NOTHROW(res->norm1(nn));
  TEUCHOS_TEST_COMPARE(nn[0], >, 1e3, out, success);

  // res contains exactly the same data as brres, the only difference is
  // that res is a MultiVector with 8 sub blocks and brres a MultiVector with
  // 3 nested sub blocks (compatible to brop)
  res->update(-STS::one(), *exp, STS::one());

  nn[0] = STS::magnitude(STS::one());
  TEST_NOTHROW(res->norm1(nn));
  TEST_EQUALITY(nn[0], STS::zero());
  TEST_NOTHROW(res->norm2(nn));
  TEST_EQUALITY(nn[0], STS::zero());
  TEST_NOTHROW(res->normInf(nn));
  TEST_EQUALITY(nn[0], STS::zero());

  // compatibility with plain maps
  Teuchos::RCP<const Map> map         = bop->getRangeMap();
  Teuchos::RCP<const BlockedMap> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap>(map);
  Teuchos::RCP<const Map> rgMap       = bmap->getFullMap();

  Teuchos::RCP<VectorClass> vrnd = VectorFactoryClass::Build(rgMap, true);
  Teuchos::RCP<VectorClass> vexp = VectorFactoryClass::Build(rgMap, true);
  Teuchos::RCP<VectorClass> vres = VectorFactoryClass::Build(rgMap, true);
  vrnd->randomize();

  // apply with plain blocked operators works with plain vectors
  TEST_NOTHROW(bop->apply(*vrnd, *vexp));

  // nested blocked operators do not work with plain vectors
  TEST_NOTHROW(brop->apply(*vrnd, *vres));

  vres->update(-STS::one(), *vexp, STS::one());
  TEUCHOS_TEST_COMPARE(vres->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(vres->normInf(), <, 5e-14, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperatorApply2, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 6 3 2 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(brop->getFullRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(brop->getFullRangeMap(), true);
  ones->putScalar(STS::one());

  brop->apply(*ones, *res);

  TEST_EQUALITY(res->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(1230)));

  // build gloabl vector with one entries (blocked version)
  Teuchos::RCP<VectorClass> bones = VectorFactoryClass::Build(brop->getDomainMap(), true);
  Teuchos::RCP<VectorClass> bres  = VectorFactoryClass::Build(brop->getRangeMap(), true);
  bones->putScalar(STS::one());

  brop->apply(*bones, *bres);

  TEST_EQUALITY(bres->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(1230)));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperatorApplyThyra, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVectorClass;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 [ 2 3 4 ] 5 ] [6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // build gloabl vector with one entries
  // The MultiVector objects "ones" and "exp" are BlockedMultiVectors with 8 sub blocks
  // compatible to bop
  Teuchos::RCP<MultiVectorClass> ones = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  Teuchos::RCP<MultiVectorClass> exp  = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  ones->putScalar(STS::one());
  bop->apply(*ones, *exp);

  // calculate 1-norm of result vector
  Teuchos::Array<typename STS::magnitudeType> nn(exp->getNumVectors());
  TEST_NOTHROW(exp->norm1(nn));
  TEST_EQUALITY(nn[0], comm->getSize() * 4485 * STS::magnitude(STS::one()));

  // overwrite result vector values
  exp->putScalar(STS::magnitude(STS::zero()));

  // reorganize "ones" and "res" to be BlockedMultiVectors with 3 sub blocks (nested)
  // compatible to brop
  // They use and work with the same 8 sub vectors from "ones" and "exp"
  Teuchos::RCP<const MultiVectorClass> cones = Teuchos::rcp_const_cast<const MultiVectorClass>(ones);
  Teuchos::RCP<const MultiVectorClass> brones =
      Xpetra::buildReorderedBlockedMultiVector(brm, Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(cones));
  Teuchos::RCP<MultiVectorClass> res        = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  Teuchos::RCP<const MultiVectorClass> cres = Teuchos::rcp_const_cast<const MultiVectorClass>(res);
  Teuchos::RCP<const MultiVectorClass> brcres =
      Xpetra::buildReorderedBlockedMultiVector(brm, Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(cres));
  Teuchos::RCP<MultiVectorClass> brres = Teuchos::rcp_const_cast<MultiVectorClass>(brcres);

  brop->apply(*brones, *brres);

  Teuchos::Array<typename STS::magnitudeType> nn2(brres->getNumVectors());
  TEST_NOTHROW(brres->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 4485 * STS::magnitude(STS::one()));
  TEST_NOTHROW(res->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 4485 * STS::magnitude(STS::one()));

  Teuchos::RCP<VectorClass> vones = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  Teuchos::RCP<VectorClass> vres  = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  vones->putScalar(STS::one());
  bop->apply(*vones, *vres);

  TEST_NOTHROW(vres->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 4485 * STS::magnitude(STS::one()));

  // not supported, yet. doImport for BlockedMultiVectors missing
  // fix this by checking the input vectors in the apply routine and switch to
  // ReorderedBlockedMultiVectors. Does this fix the problem? Do it in a two stage fashion.
  TEST_NOTHROW(brop->apply(*vones, *vres));
  TEST_NOTHROW(vres->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 4485 * STS::magnitude(STS::one()));

  Teuchos::RCP<VectorClass> vrnd2 = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  Teuchos::RCP<VectorClass> vexp2 = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  Teuchos::RCP<VectorClass> vres2 = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  vrnd2->randomize();

  // apply with plain blocked operators works with plain vectors
  TEST_NOTHROW(bop->apply(*vrnd2, *vexp2));

  // nested blocked operators do not work with plain vectors
  TEST_NOTHROW(brop->apply(*vrnd2, *vres2));

  vres2->update(-STS::one(), *vexp2, STS::one());
  TEUCHOS_TEST_COMPARE(vres2->norm2(), <, 5e-14, out, success);
  TEUCHOS_TEST_COMPARE(vres2->normInf(), <, 5e-14, out, success);

  // build gloabl vector with one entries (blocked version)
  Teuchos::RCP<VectorClass> bones = VectorFactoryClass::Build(brop->getDomainMap(), true);
  Teuchos::RCP<VectorClass> bres  = VectorFactoryClass::Build(brop->getRangeMap(), true);
  bones->putScalar(STS::one());

  brop->apply(*bones, *bres);

  TEST_EQUALITY(bres->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(4485)));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperatorApplyThyraSmall, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVectorClass;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 3;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 2 ] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // build gloabl vector with one entries
  // The MultiVector objects "ones" and "exp" are BlockedMultiVectors with 8 sub blocks
  // compatible to bop
  Teuchos::RCP<MultiVectorClass> ones = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  Teuchos::RCP<MultiVectorClass> exp  = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  ones->putScalar(STS::one());
  bop->apply(*ones, *exp);

  // calculate 1-norm of result vector
  Teuchos::Array<typename STS::magnitudeType> nn(exp->getNumVectors());
  TEST_NOTHROW(exp->norm1(nn));
  TEST_EQUALITY(nn[0], comm->getSize() * 45 * STS::magnitude(STS::one()));

  // overwrite result vector values
  exp->putScalar(STS::magnitude(STS::zero()));

  // exp->describe(out,Teuchos::VERB_EXTREME);

  // reorganize "ones" and "res" to be BlockedMultiVectors with 2 sub blocks (nested)
  // compatible to brop
  Teuchos::RCP<const MultiVectorClass> cones = Teuchos::rcp_const_cast<const MultiVectorClass>(ones);
  Teuchos::RCP<const MultiVectorClass> brones =
      Xpetra::buildReorderedBlockedMultiVector(brm, Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(cones));
  Teuchos::RCP<MultiVectorClass> res        = MultiVectorFactoryClass::Build(bop->getRangeMap(), 1, true);
  Teuchos::RCP<const MultiVectorClass> cres = Teuchos::rcp_const_cast<const MultiVectorClass>(res);
  Teuchos::RCP<const MultiVectorClass> brcres =
      Xpetra::buildReorderedBlockedMultiVector(brm, Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(cres));
  Teuchos::RCP<MultiVectorClass> brres = Teuchos::rcp_const_cast<MultiVectorClass>(brcres);

  // Note: the result is both contained in res and brres.
  //       brres uses the same underlying vectors!
  brop->apply(*brones, *brres);

  Teuchos::Array<typename STS::magnitudeType> nn2(brres->getNumVectors());
  TEST_NOTHROW(brres->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 45 * STS::magnitude(STS::one()));
  TEST_NOTHROW(res->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 45 * STS::magnitude(STS::one()));

  Teuchos::RCP<VectorClass> vones = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  Teuchos::RCP<VectorClass> vres  = VectorFactoryClass::Build(bop->getFullRangeMap(), true);
  vones->putScalar(STS::one());
  bop->apply(*vones, *vres);

  TEST_NOTHROW(vres->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 45 * STS::magnitude(STS::one()));

  // not supported, yet. doImport for BlockedMultiVectors missing
  // fix this by checking the input vectors in the apply routine and switch to
  // ReorderedBlockedMultiVectors. Does this fix the problem? Do it in a two stage fashion.
  TEST_NOTHROW(brop->apply(*vones, *vres));

  TEST_NOTHROW(vres->norm1(nn2));
  TEST_EQUALITY(nn2[0], comm->getSize() * 45 * STS::magnitude(STS::one()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReorderBlockOperatorApply2Thyra, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 8;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 6 3 2 ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(brop->getFullRangeMap(), true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(brop->getFullRangeMap(), true);
  ones->putScalar(STS::one());

  brop->apply(*ones, *res);

  TEST_EQUALITY(res->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(1230)));

  // build gloabl vector with one entries (blocked version)
  Teuchos::RCP<VectorClass> bones = VectorFactoryClass::Build(brop->getDomainMap(), true);
  Teuchos::RCP<VectorClass> bres  = VectorFactoryClass::Build(brop->getRangeMap(), true);
  bones->putScalar(STS::one());

  brop->apply(*bones, *bres);

  TEST_EQUALITY(bres->norm1(), STS::magnitude(Teuchos::as<Scalar>(comm->getSize()) * Teuchos::as<Scalar>(1230)));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ReadWriteBlockedMatrix, M, MA, Scalar, LO, GO, Node) {
  // TODO: it seems that the Tpetra matrix reader is only working for standard maps??

  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar, LO, GO, Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar, LO, GO, Node> MapExtractorFactoryClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  // generate problem
  GO nEle                                = 63;
  const Teuchos::RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

  LO NumMyElements                              = map->getLocalNumElements();
  GO NumGlobalElements                          = map->getGlobalNumElements();
  Teuchos::ArrayView<const GO> MyGlobalElements = map->getLocalElementList();

  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node>> A =
      Xpetra::CrsMatrixFactory<Scalar, LO, GO, Node>::Build(map, 3);
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == true || A->isFillActive() == false, std::runtime_error, "");

  for (LO i = 0; i < NumMyElements; i++) {
    if (MyGlobalElements[i] == 0) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(Teuchos::as<Scalar>(i) * STS::one(), -1.0));
    } else if (MyGlobalElements[i] == NumGlobalElements - 1) {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i]),
                            Teuchos::tuple<Scalar>(-1.0, Teuchos::as<Scalar>(i) * STS::one()));
    } else {
      A->insertGlobalValues(MyGlobalElements[i],
                            Teuchos::tuple<GO>(MyGlobalElements[i] - 1, MyGlobalElements[i], MyGlobalElements[i] + 1),
                            Teuchos::tuple<Scalar>(-1.0, Teuchos::as<Scalar>(i) * STS::one(), -1.0));
    }
  }

  A->fillComplete();
  TEUCHOS_TEST_FOR_EXCEPTION(A->isFillComplete() == false || A->isFillActive() == true, std::runtime_error, "");

  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> mat =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>(A));

  Teuchos::Array<GO> gids1;
  Teuchos::Array<GO> gids2;
  for (LO i = 0; i < NumMyElements; i++) {
    if (i % 3 < 2)
      gids1.push_back(map->getGlobalElement(i));
    else
      gids2.push_back(map->getGlobalElement(i));
  }

  const Teuchos::RCP<const MapClass> map1 = MapFactoryClass::Build(lib,
                                                                   Teuchos::OrdinalTraits<GO>::invalid(),
                                                                   gids1.view(0, gids1.size()),
                                                                   0,
                                                                   comm);
  const Teuchos::RCP<const MapClass> map2 = MapFactoryClass::Build(lib,
                                                                   Teuchos::OrdinalTraits<GO>::invalid(),
                                                                   gids2.view(0, gids2.size()),
                                                                   0,
                                                                   comm);

  std::vector<Teuchos::RCP<const MapClass>> xmaps;
  xmaps.push_back(map1);
  xmaps.push_back(map2);

  Teuchos::RCP<const MapExtractorClass> rowMapExtractormap_extractor = MapExtractorFactoryClass::Build(map, xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bMat =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*mat, rowMapExtractormap_extractor, rowMapExtractormap_extractor);

  // Write matrices out, read fine A back in, and check that the read was ok
  // by using a matvec with a random vector.
  // JJH: 22-Feb-2016 Append scalar type to file name. The theory is that for dashboard
  //      tests with multiple Scalar instantiations of this test, a test with Scalar type
  //      A could try to read in the results of the test with Scalar type B, simply because
  //      the test with type B overwrote A's output matrix file.  A better solution would be
  //      to write to a file stream, but this would involve writing new interfaces to Epetra's
  //      file I/O capabilities.
  std::string tname = "BLOCKEDMATRIX";
  tname             = tname + typeid(Scalar).name();
  tname             = tname + typeid(LO).name();
  tname             = tname + typeid(GO).name();
#ifdef HAVE_MUELU_KOKKOS
  std::string nn = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<typename Node::execution_space>::name();
  nn.erase(std::remove(nn.begin(), nn.end(), '/'), nn.end());
  tname = tname + nn;
#endif
  tname = "_" + tname;

  const bool writeAllMaps = true;
  Xpetra::IO<Scalar, LO, GO, Node>::WriteBlockedCrsMatrix(tname, *bMat, writeAllMaps);
  Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bMat2 = Xpetra::IO<Scalar, LO, GO, Node>::ReadBlockedCrsMatrix(tname, lib, comm);

  TEST_EQUALITY(bMat->getMatrix(0, 0)->getGlobalNumEntries(), bMat2->getMatrix(0, 0)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(0, 1)->getGlobalNumEntries(), bMat2->getMatrix(0, 1)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1, 0)->getGlobalNumEntries(), bMat2->getMatrix(1, 0)->getGlobalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1, 1)->getGlobalNumEntries(), bMat2->getMatrix(1, 1)->getGlobalNumEntries());

  TEST_EQUALITY(bMat->getMatrix(0, 0)->getLocalNumEntries(), bMat2->getMatrix(0, 0)->getLocalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(0, 1)->getLocalNumEntries(), bMat2->getMatrix(0, 1)->getLocalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1, 0)->getLocalNumEntries(), bMat2->getMatrix(1, 0)->getLocalNumEntries());
  TEST_EQUALITY(bMat->getMatrix(1, 1)->getLocalNumEntries(), bMat2->getMatrix(1, 1)->getLocalNumEntries());

  TEST_EQUALITY(bMat->getMatrix(0, 0)->getFrobeniusNorm(), bMat2->getMatrix(0, 0)->getFrobeniusNorm());
  TEST_EQUALITY(bMat->getMatrix(0, 1)->getFrobeniusNorm(), bMat2->getMatrix(0, 1)->getFrobeniusNorm());
  TEST_EQUALITY(bMat->getMatrix(1, 0)->getFrobeniusNorm(), bMat2->getMatrix(1, 0)->getFrobeniusNorm());
  TEST_EQUALITY(bMat->getMatrix(1, 1)->getFrobeniusNorm(), bMat2->getMatrix(1, 1)->getFrobeniusNorm());

  TEST_EQUALITY(bMat->getRangeMapExtractor()->getMap(0)->isSameAs(*(bMat2->getRangeMapExtractor()->getMap(0))), true);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getMap(0)->isSameAs(*(bMat2->getDomainMapExtractor()->getMap(0))), true);

  TEST_EQUALITY(bMat->getRangeMapExtractor()->getFullMap()->isSameAs(*(bMat2->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getFullMap()->isSameAs(*(bMat2->getDomainMapExtractor()->getFullMap())), true);

  // these tests are false with Tpetra? TODO check me: why only in Tpetra?
  // bMat2 is always in Xpetra mode so far. This is, since the Read routine and Write routine for the MapExtractor do not really
  // consider the Thyra mode so far.
  // TEST_EQUALITY(bMat->getRangeMapExtractor()->getMap(1)->isSameAs(*(bMat2->getRangeMapExtractor()->getMap(1))),true);
  // TEST_EQUALITY(bMat->getDomainMapExtractor()->getMap(1)->isSameAs(*(bMat2->getDomainMapExtractor()->getMap(1))),true);

  TEST_EQUALITY(bMat->getMatrix(0, 0)->getRowMap()->isSameAs(*(bMat2->getMatrix(0, 0)->getRowMap())), true);
  TEST_EQUALITY(bMat->getMatrix(0, 1)->getRowMap()->isSameAs(*(bMat2->getMatrix(0, 1)->getRowMap())), true);
  TEST_EQUALITY(bMat->getMatrix(1, 0)->getRowMap()->isSameAs(*(bMat2->getMatrix(1, 0)->getRowMap())), true);
  TEST_EQUALITY(bMat->getMatrix(1, 1)->getRowMap()->isSameAs(*(bMat2->getMatrix(1, 1)->getRowMap())), true);

  TEST_EQUALITY(bMat->getMatrix(0, 0)->getColMap()->isSameAs(*(bMat2->getMatrix(0, 0)->getColMap())), true);
  TEST_EQUALITY(bMat->getMatrix(0, 1)->getColMap()->isSameAs(*(bMat2->getMatrix(0, 1)->getColMap())), true);
  // the following test fails with Teptra. Why?
  // TEST_EQUALITY(bMat->getMatrix(1,0)->getColMap()->isSameAs(*(bMat2->getMatrix(1,0)->getColMap())),true);
  TEST_EQUALITY(bMat->getMatrix(1, 1)->getColMap()->isSameAs(*(bMat2->getMatrix(1, 1)->getColMap())), true);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones_A   = VectorFactoryClass::Build(bMat->getRangeMap(), true);
  Teuchos::RCP<VectorClass> exp      = VectorFactoryClass::Build(bMat->getRangeMap(), true);
  Teuchos::RCP<VectorClass> ones_bOp = VectorFactoryClass::Build(bMat2->getRangeMap(), true);
  Teuchos::RCP<VectorClass> res      = VectorFactoryClass::Build(bMat2->getRangeMap(), true);
  ones_A->putScalar(STS::one());
  ones_bOp->putScalar(STS::one());

  bMat->apply(*ones_A, *exp);
  bMat2->apply(*ones_bOp, *res);

  TEST_EQUALITY(res->norm2(), exp->norm2());
  TEST_EQUALITY(res->normInf(), exp->normInf());
  TEST_EQUALITY(bMat->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat->getDomainMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat2->getRangeMapExtractor()->getThyraMode(), false);  // thyra mode is not correctly transferred!!
  TEST_EQUALITY(bMat2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bMat2->getRangeMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat2->getRangeMap(1)->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(bMat2->getDomainMap(0)->getMinAllGlobalIndex(), 0);
  TEST_INEQUALITY(bMat2->getDomainMap(1)->getMinAllGlobalIndex(), 0);
}

/// simple test routine for the apply function of BlockedCrsMatrix
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, Apply, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::StridedMap<LO, GO, Node> StridedMapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::StridedMapFactory<LO, GO, Node> StridedMapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar, LO, GO, Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar, LO, GO, Node> MapExtractorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const Teuchos::RCP<const MapClass> pointmap = MapFactoryClass::Build(lib, 12, 0, comm);

  // generate local maps for loading matrices
  Teuchos::Array<GO> velgidvec;  // global strided maps
  Teuchos::Array<GO> pregidvec;
  Teuchos::Array<GO> fullgidvec;  // full global map
  for (LO i = 0; i < Teuchos::as<LO>(pointmap->getLocalNumElements()); i++) {
    // loop over all local ids in pointmap

    // get corresponding global id
    GO gid = pointmap->getGlobalElement(i);

    // store global strided gids
    velgidvec.push_back(3 * gid);
    velgidvec.push_back(3 * gid + 1);
    pregidvec.push_back(3 * gid + 2);

    // gid for full map
    fullgidvec.push_back(3 * gid);
    fullgidvec.push_back(3 * gid + 1);
    fullgidvec.push_back(3 * gid + 2);
  }

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  const Teuchos::RCP<const StridedMapClass> velmap  = StridedMapFactoryClass::Build(lib,
                                                                                    Teuchos::OrdinalTraits<GO>::invalid(),
                                                                                    velgidvec(),
                                                                                    0,
                                                                                    stridingInfo,
                                                                                    comm, 0);
  const Teuchos::RCP<const StridedMapClass> premap  = StridedMapFactoryClass::Build(lib,
                                                                                    Teuchos::OrdinalTraits<GO>::invalid(),
                                                                                    pregidvec(),
                                                                                    0,
                                                                                    stridingInfo,
                                                                                    comm, 1);
  const Teuchos::RCP<const StridedMapClass> fullmap = StridedMapFactoryClass::Build(lib,
                                                                                    Teuchos::OrdinalTraits<GO>::invalid(),
                                                                                    fullgidvec(),
                                                                                    0,
                                                                                    stridingInfo,
                                                                                    comm, -1);

  std::string tname = typeid(Scalar).name();
  if (tname.find("complex") != std::string::npos) {
    std::cout << "Skip test for scalar=" << tname << std::endl;
    return;
  }

  Teuchos::RCP<MatrixClass> A = Xpetra::IO<Scalar, LO, GO, Node>::Read("A.mat", fullmap->getMap());

  std::vector<Teuchos::RCP<const MapClass>> xmaps;
  xmaps.push_back(velmap);
  xmaps.push_back(premap);

  Teuchos::RCP<const MapExtractorClass> rowMapExtractormap_extractor = MapExtractorFactoryClass::Build(fullmap, xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bOp =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*A, rowMapExtractormap_extractor, rowMapExtractormap_extractor);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(fullmap, true);
  ones->putScalar(STS::one());
  rnd->randomize();

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  magnitudeType tol1  = Teuchos::ScalarTraits<magnitudeType>::eps();
  magnitudeType tol2  = 500 * tol1;

  A->apply(*ones, *exp);
  bOp->apply(*ones, *res);
  res->update(-STS::one(), *exp, STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, tol1, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, tol1, out, success);

  A->apply(*rnd, *exp);
  bOp->apply(*rnd, *res);
  res->update(-STS::one(), *exp, STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, tol2, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, tol2, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, getLocalDiagCopy, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMapClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVectorClass;
  typedef Xpetra::BlockedVector<Scalar, LO, GO, Node> BlockedVectorClass;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 5;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<VectorClass> vorig = VectorFactoryClass::Build(bop->getRangeMap(), true);

  bop->getLocalDiagCopy(*vorig);

  Teuchos::RCP<BlockedVectorClass> bvorig = Teuchos::rcp_dynamic_cast<BlockedVectorClass>(vorig);
  TEST_EQUALITY(bvorig.is_null(), false);
  TEST_EQUALITY(bvorig->getBlockedMap()->getNumMaps(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const MultiVectorClass> mvorig = bvorig->Merge();
  TEST_EQUALITY(mvorig->getMap()->getMaxAllGlobalIndex(), bvorig->getMap()->getMaxAllGlobalIndex());
  TEST_EQUALITY(mvorig->getNumVectors(), 1);

  Teuchos::ArrayRCP<const Scalar> vdataorig = mvorig->getData(0);
  bool bCheck                               = true;
  for (int i = 0; i < 5; i++)
    if (vdataorig[i] != Teuchos::as<Scalar>(1.0)) bCheck = false;
  for (int i = 5; i < 10; i++)
    if (vdataorig[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for (int i = 10; i < 20; i++)
    if (vdataorig[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  for (int i = 20; i < 40; i++)
    if (vdataorig[i] != Teuchos::as<Scalar>(4.0)) bCheck = false;
  for (int i = 40; i < 80; i++)
    if (vdataorig[i] != Teuchos::as<Scalar>(5.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // reordered blocked operator (Xpetra style)
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 4 [3 2] 1 0]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  GO goNumRows = Teuchos::as<GO>(Teuchos::ScalarTraits<GO>::pow(2, noBlocks - 2)) * 10 * comm->getSize();

  TEST_EQUALITY(brop->Rows(), 4);
  TEST_EQUALITY(brop->Cols(), 4);
  TEST_EQUALITY(brop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(brop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));

  Teuchos::RCP<const MapClass> map         = brop->getRangeMap();
  Teuchos::RCP<const BlockedMapClass> bmap = Teuchos::rcp_dynamic_cast<const BlockedMapClass>(map);
  TEST_EQUALITY(bmap.is_null(), false);
  TEST_EQUALITY(bmap->getNumMaps(), 4);
  Teuchos::RCP<VectorClass> v         = VectorFactoryClass::Build(brop->getRangeMap(), true);
  Teuchos::RCP<BlockedVectorClass> bv = Teuchos::rcp_dynamic_cast<BlockedVectorClass>(v);
  TEST_EQUALITY(bv.is_null(), false);

  brop->getLocalDiagCopy(*v);

  mvorig = bv->Merge();
  TEST_EQUALITY(mvorig->getMap()->getMaxAllGlobalIndex(), v->getMap()->getMaxAllGlobalIndex());
  TEST_EQUALITY(mvorig->getNumVectors(), 1);

  Teuchos::ArrayRCP<const Scalar> vdata = mvorig->getData(0);
  bCheck                                = true;
  for (int i = 0; i < 40; i++)
    if (vdata[i] != Teuchos::as<Scalar>(5.0)) bCheck = false;
  for (int i = 40; i < 60; i++)
    if (vdata[i] != Teuchos::as<Scalar>(4.0)) bCheck = false;
  for (int i = 60; i < 70; i++)
    if (vdata[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  for (int i = 70; i < 75; i++)
    if (vdata[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for (int i = 75; i < 80; i++)
    if (vdata[i] != Teuchos::as<Scalar>(1.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // Thyra style (reordered) operator
  Teuchos::RCP<const BlockedCrsMatrixClass> btop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brtop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, btop));

  TEST_EQUALITY(brtop->Rows(), 4);
  TEST_EQUALITY(brtop->Cols(), 4);
  TEST_EQUALITY(brtop->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));
  TEST_EQUALITY(brtop->getDomainMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(goNumRows));

  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(brtop->getRangeMap(), true);

  brtop->getLocalDiagCopy(*v2);

  Teuchos::RCP<BlockedVectorClass> bv2 = Teuchos::rcp_dynamic_cast<BlockedVectorClass>(v2);
  TEST_EQUALITY(bv2.is_null(), false);

  bCheck             = true;
  int expectedResult = 5;
  for (int k = 0; k < 4; k++) {
    Teuchos::RCP<MultiVectorClass> subvec = bv2->getMultiVector(k);
    TEST_EQUALITY(subvec.is_null(), false);
    TEST_EQUALITY(subvec->getNumVectors(), 1);
    TEST_EQUALITY(Teuchos::rcp_dynamic_cast<VectorClass>(subvec) != Teuchos::null, true);
    Teuchos::RCP<BlockedMultiVectorClass> bsubvec = Teuchos::rcp_dynamic_cast<BlockedMultiVectorClass>(subvec);
    if (bsubvec == Teuchos::null) {
      Teuchos::ArrayRCP<const Scalar> vdata2 = subvec->getData(0);
      for (size_t l = 0; l < Teuchos::as<size_t>(vdata2.size()); l++) {
        if (vdata2[l] != Teuchos::as<Scalar>(expectedResult)) bCheck = false;
      }
      expectedResult--;
    } else {
      for (size_t m = 0; m < bsubvec->getBlockedMap()->getNumMaps(); m++) {
        Teuchos::RCP<MultiVectorClass> ssubvec = bsubvec->getMultiVector(m);
        TEST_EQUALITY(ssubvec->getNumVectors(), 1);
        TEST_EQUALITY(Teuchos::rcp_dynamic_cast<BlockedMultiVectorClass>(ssubvec) == Teuchos::null, true);
        Teuchos::ArrayRCP<const Scalar> vdata3 = ssubvec->getData(0);
        for (size_t l2 = 0; l2 < Teuchos::as<size_t>(vdata3.size()); l2++) {
          if (vdata3[l2] != Teuchos::as<Scalar>(expectedResult)) bCheck = false;
        }
        expectedResult--;
      }
    }
  }
  TEST_EQUALITY(bCheck, true);

  mvorig                                 = bv2->Merge();
  Teuchos::ArrayRCP<const Scalar> vdata2 = mvorig->getData(0);
  bCheck                                 = true;
  for (int i = 0; i < 40; i++)
    if (vdata2[i] != Teuchos::as<Scalar>(5.0)) bCheck = false;
  for (int i = 40; i < 60; i++)
    if (vdata2[i] != Teuchos::as<Scalar>(4.0)) bCheck = false;
  for (int i = 60; i < 70; i++)
    if (vdata2[i] != Teuchos::as<Scalar>(3.0)) bCheck = false;
  for (int i = 70; i < 75; i++)
    if (vdata2[i] != Teuchos::as<Scalar>(2.0)) bCheck = false;
  for (int i = 75; i < 80; i++)
    if (vdata2[i] != Teuchos::as<Scalar>(1.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, leftScale, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 3;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  Teuchos::RCP<VectorClass> v1 = VectorFactoryClass::Build(bop->getRangeMap(), true);
  bop->getLocalDiagCopy(*v1);

  Teuchos::RCP<BlockedCrsMatrixClass> bop_nonconst = Teuchos::rcp_const_cast<BlockedCrsMatrixClass>(bop);

  Teuchos::RCP<VectorClass> s = VectorFactoryClass::Build(bop_nonconst->getRangeMap(), true);
  s->putScalar(Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(2.0));

  bop_nonconst->leftScale(*s);

  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(bop_nonconst->getRangeMap(), true);
  bop_nonconst->getLocalDiagCopy(*v2);

  TEST_EQUALITY(v1->norm1() * Teuchos::as<Scalar>(2.0), v2->norm1());
  TEST_EQUALITY(v1->norm2() * Teuchos::as<Scalar>(2.0), v2->norm2());
  TEST_EQUALITY(v1->normInf() * Teuchos::as<Scalar>(2.0), v2->normInf());
  v1->update(Teuchos::as<Scalar>(-0.5), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEST_EQUALITY(v1->norm1(), Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()));

  // reordered Xpetra operator
  bop                                                 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(8, *comm);
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 [ 2 3 4 ] 5 ] [6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  Teuchos::RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop_nonconst =
      Teuchos::rcp_const_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop);

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);

  s = VectorFactoryClass::Build(brop->getRangeMap(), true);
  s->putScalar(Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(2.0));

  v1 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v1);

  brop_nonconst->leftScale(*s);

  v2 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v2);

  TEST_EQUALITY(v1->norm1() * Teuchos::as<Scalar>(2.0), v2->norm1());
  TEST_EQUALITY(v1->norm2() * Teuchos::as<Scalar>(2.0), v2->norm2());
  TEST_EQUALITY(v1->normInf() * Teuchos::as<Scalar>(2.0), v2->normInf());
  v1->update(Teuchos::as<Scalar>(-0.5), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEST_EQUALITY(v1->norm1(), Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()));

  // reordered Thyra operator
  bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(8, *comm);

  brop          = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  brop_nonconst = Teuchos::rcp_const_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop);

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  s = VectorFactoryClass::Build(brop->getRangeMap(), true);
  s->putScalar(Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(2.0));

  v1 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v1);

  brop_nonconst->leftScale(*s);

  v2 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v2);

  TEST_EQUALITY(v1->norm1() * Teuchos::as<Scalar>(2.0), v2->norm1());
  TEST_EQUALITY(v1->norm2() * Teuchos::as<Scalar>(2.0), v2->norm2());
  TEST_EQUALITY(v1->normInf() * Teuchos::as<Scalar>(2.0), v2->normInf());
  v1->update(Teuchos::as<Scalar>(-0.5), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEST_EQUALITY(v1->norm1(), Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, rightScale, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 3;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<VectorClass> v1 = VectorFactoryClass::Build(bop->getRangeMap(), true);
  bop->getLocalDiagCopy(*v1);

  Teuchos::RCP<BlockedCrsMatrixClass> bop_nonconst = Teuchos::rcp_const_cast<BlockedCrsMatrixClass>(bop);

  Teuchos::RCP<VectorClass> s = VectorFactoryClass::Build(bop_nonconst->getDomainMap(), true);
  s->putScalar(Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(2.0));

  bop_nonconst->rightScale(*s);

  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(bop_nonconst->getRangeMap(), true);
  bop_nonconst->getLocalDiagCopy(*v2);

  /*Teuchos::ArrayRCP< const Scalar > v1d = v1->getData(0);
  Teuchos::ArrayRCP< const Scalar > v2d = v2->getData(0);
  bool bCheck = true;
  for(int i=0; i<20; i++)  if(v1d[i] * Teuchos::as<Scalar>(2.0) != v2d[i]) bCheck = false;
  TEST_EQUALITY(bCheck, true);*/
  TEST_EQUALITY(v1->norm1() * Teuchos::as<Scalar>(2.0), v2->norm1());
  TEST_EQUALITY(v1->norm2() * Teuchos::as<Scalar>(2.0), v2->norm2());
  TEST_EQUALITY(v1->normInf() * Teuchos::as<Scalar>(2.0), v2->normInf());
  v1->update(Teuchos::as<Scalar>(-0.5), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEST_EQUALITY(v1->norm1(), Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()));

  // reordered Xpetra operator
  bop                                                 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(8, *comm);
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [ 1 [ 2 3 4 ] 5 ] [6 7] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  Teuchos::RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop_nonconst =
      Teuchos::rcp_const_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop);

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);

  s = VectorFactoryClass::Build(brop->getDomainMap(), true);
  s->putScalar(Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(2.0));

  v1 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v1);

  brop_nonconst->rightScale(*s);

  v2 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v2);

  /*v1d = v1->getData(0);
  v2d = v2->getData(0);
  bCheck = true;
  for(int i=0; i<640; i++)  if(v1d[i] * Teuchos::as<Scalar>(2.0) != v2d[i]) bCheck = false;
  TEST_EQUALITY(bCheck, true);*/
  TEST_EQUALITY(v1->norm1() * Teuchos::as<Scalar>(2.0), v2->norm1());
  TEST_EQUALITY(v1->norm2() * Teuchos::as<Scalar>(2.0), v2->norm2());
  TEST_EQUALITY(v1->normInf() * Teuchos::as<Scalar>(2.0), v2->normInf());
  v1->update(Teuchos::as<Scalar>(-0.5), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEST_EQUALITY(v1->norm1(), Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()));

  // reordered Thyra operator
  bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(8, *comm);

  brop          = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  brop_nonconst = Teuchos::rcp_const_cast<Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(brop);

  TEST_EQUALITY(brop->Rows(), 3);
  TEST_EQUALITY(brop->Cols(), 3);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);

  s = VectorFactoryClass::Build(brop->getDomainMap(), true);
  s->putScalar(Teuchos::ScalarTraits<Scalar>::one() * Teuchos::as<Scalar>(2.0));

  v1 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v1);

  brop_nonconst->rightScale(*s);

  v2 = VectorFactoryClass::Build(brop_nonconst->getRangeMap(), true);
  brop_nonconst->getLocalDiagCopy(*v2);

  /*v1d = v1->getData(0);
  v2d = v2->getData(0);
  bCheck = true;
  for(int i=0; i<640; i++)  if(v1d[i] * Teuchos::as<Scalar>(2.0) != v2d[i]) bCheck = false;
  TEST_EQUALITY(bCheck, true);*/
  TEST_EQUALITY(v1->norm1() * Teuchos::as<Scalar>(2.0), v2->norm1());
  TEST_EQUALITY(v1->norm2() * Teuchos::as<Scalar>(2.0), v2->norm2());
  TEST_EQUALITY(v1->normInf() * Teuchos::as<Scalar>(2.0), v2->normInf());
  v1->update(Teuchos::as<Scalar>(-0.5), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEST_EQUALITY(v1->norm1(), Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()));
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, deepCopy, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::MatrixFactory<Scalar, LO, GO, Node> MatrixFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 5;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const MatrixClass> A2             = MatrixFactoryClass::BuildCopy(bop);
  Teuchos::RCP<const BlockedCrsMatrixClass> bop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A2);
  TEST_EQUALITY(bop2->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop2->Cols(), Teuchos::as<size_t>(noBlocks));

  TEST_EQUALITY(bop2->getRangeMapExtractor()->NumMaps(), bop->getRangeMapExtractor()->NumMaps());
  TEST_EQUALITY(bop2->getDomainMapExtractor()->NumMaps(), bop->getDomainMapExtractor()->NumMaps());
  TEST_EQUALITY(bop2->getGlobalMaxNumRowEntries(), bop->getGlobalMaxNumRowEntries());
  TEST_EQUALITY(bop2->getGlobalNumEntries(), bop->getGlobalNumEntries());
  TEST_EQUALITY(bop2->getGlobalNumRows(), bop->getGlobalNumRows());
  TEST_EQUALITY(bop2->getGlobalNumCols(), bop->getGlobalNumCols());

  Teuchos::RCP<VectorClass> v1 = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(bop2->getRangeMap(), true);
  bop->getLocalDiagCopy(*v1);
  bop2->getLocalDiagCopy(*v2);

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  magnitudeType tol   = Teuchos::ScalarTraits<magnitudeType>::eps();

  v1->update(-Teuchos::ScalarTraits<Scalar>::one(), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEUCHOS_TEST_COMPARE(v1->norm2(), <, tol, out, success);
  TEUCHOS_TEST_COMPARE(v1->normInf(), <, tol, out, success);

  v1 = Teuchos::null;
  v2 = Teuchos::null;

  bop = Teuchos::null;

  TEST_EQUALITY(bop2->getRangeMapExtractor()->NumMaps(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop2->getDomainMapExtractor()->NumMaps(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop2->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bop2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bop2->getGlobalMaxNumRowEntries(), 1);
  TEST_EQUALITY(bop2->getGlobalNumRows(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 80));
  TEST_EQUALITY(bop2->getGlobalNumCols(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 80));
  TEST_EQUALITY(bop2->getMatrix(0, 0) != Teuchos::null, true);
  TEST_EQUALITY(bop2->getMatrix(0, 0)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));

  // Thyra blocked operator
  Teuchos::RCP<const BlockedCrsMatrixClass> bop3 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop3->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop3->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const MatrixClass> A4             = MatrixFactoryClass::BuildCopy(bop3);
  Teuchos::RCP<const BlockedCrsMatrixClass> bop4 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A4);
  TEST_EQUALITY(bop4->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop4->Cols(), Teuchos::as<size_t>(noBlocks));

  bop3 = Teuchos::null;

  TEST_EQUALITY(bop4->getRangeMapExtractor()->NumMaps(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop4->getDomainMapExtractor()->NumMaps(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop4->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bop4->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bop4->getGlobalMaxNumRowEntries(), 1);
  TEST_EQUALITY(bop4->getGlobalNumRows(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 80));
  TEST_EQUALITY(bop4->getGlobalNumCols(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 80));
  TEST_EQUALITY(bop4->getMatrix(0, 0) != Teuchos::null, true);
  TEST_EQUALITY(bop4->getMatrix(0, 0)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 5));

  // Nested Xpetra blocked operator
  bop                                                 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 1 [ [ 2 4 0 ] 3] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);

  Teuchos::RCP<const MatrixClass> A               = MatrixFactoryClass::BuildCopy(brop);
  Teuchos::RCP<const BlockedCrsMatrixClass> brop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A);
  brop                                            = Teuchos::null;
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  Teuchos::RCP<const BlockedCrsMatrixClass> brop200   = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(0, 0));
  Teuchos::RCP<const BlockedCrsMatrixClass> brop211   = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(1, 1));
  Teuchos::RCP<const BlockedCrsMatrixClass> brop21100 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop211->getMatrix(0, 0));
  TEST_EQUALITY(brop200->Rows(), 1);
  TEST_EQUALITY(brop200->Cols(), 1);
  TEST_EQUALITY(brop211->Rows(), 2);
  TEST_EQUALITY(brop211->Cols(), 2);
  TEST_EQUALITY(brop21100->Rows(), 3);
  TEST_EQUALITY(brop21100->Cols(), 3);
  TEST_EQUALITY(brop21100->getMatrix(0, 0)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(brop21100->getMatrix(1, 1)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop200->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop200->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop211->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop211->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop21100->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop21100->getDomainMapExtractor()->getThyraMode(), false);

  // Nested Thyra blocked operator
  bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  brop = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);

  A     = MatrixFactoryClass::BuildCopy(brop);
  brop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A);
  brop  = Teuchos::null;
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  brop200   = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(0, 0));
  brop211   = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop2->getMatrix(1, 1));
  brop21100 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(brop211->getMatrix(0, 0));
  TEST_EQUALITY(brop200->Rows(), 1);
  TEST_EQUALITY(brop200->Cols(), 1);
  TEST_EQUALITY(brop211->Rows(), 2);
  TEST_EQUALITY(brop211->Cols(), 2);
  TEST_EQUALITY(brop21100->Rows(), 3);
  TEST_EQUALITY(brop21100->Cols(), 3);
  TEST_EQUALITY(brop21100->getMatrix(0, 0)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 10));
  TEST_EQUALITY(brop21100->getMatrix(1, 1)->getRangeMap()->getGlobalNumElements(), Teuchos::as<Xpetra::global_size_t>(comm->getSize() * 40));
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop200->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop200->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop211->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop211->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop21100->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop21100->getDomainMapExtractor()->getThyraMode(), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, transformThyra2XpetraGIDs, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::MapUtils<LO, GO, Node> MapUtilsClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  Teuchos::Array<GO> ovltGIDs;
  Teuchos::Array<GO> novltGIDs;
  Teuchos::Array<GO> novlxGIDs;

  for (int i = 0; i < 10; i++) {
    novltGIDs.append(comm->getRank() * 10 + Teuchos::as<GO>(i));
    novlxGIDs.append(comm->getRank() * 40 + Teuchos::as<GO>(i) * 10 + 111);
    ovltGIDs.append(comm->getRank() * 10 + Teuchos::as<GO>(i));
  }
  if (comm->getRank() > 0) ovltGIDs.append(comm->getRank() * 10 - 1);
  if (comm->getRank() < comm->getSize() - 1) ovltGIDs.append(comm->getRank() * 10 + 10);

  Teuchos::RCP<const MapClass> ovltMap  = MapFactoryClass::Build(lib, Teuchos::OrdinalTraits<GO>::invalid(), ovltGIDs(), 0, comm);
  Teuchos::RCP<const MapClass> novltMap = MapFactoryClass::Build(lib, Teuchos::OrdinalTraits<GO>::invalid(), novltGIDs(), 0, comm);
  Teuchos::RCP<const MapClass> novlxMap = MapFactoryClass::Build(lib, Teuchos::OrdinalTraits<GO>::invalid(), novlxGIDs(), 0, comm);

  Teuchos::RCP<MapClass> xmap = MapUtilsClass::transformThyra2XpetraGIDs(
      *ovltMap,
      *novltMap,
      *novlxMap);

  TEST_EQUALITY(xmap->getLocalNumElements(), ovltMap->getLocalNumElements());
  TEST_EQUALITY(xmap->getGlobalNumElements(), ovltMap->getGlobalNumElements());
  for (int i = 0; i < 10; i++) {
    GO gid = xmap->getGlobalElement(i);
    TEST_EQUALITY(gid, comm->getRank() * 40 + Teuchos::as<GO>(i) * 10 + 111);
  }
  if (comm->getRank() > 0 && comm->getRank() < comm->getSize() - 2)
    TEST_EQUALITY(xmap->getGlobalElement(10), (comm->getRank() - 1) * 40 + Teuchos::as<GO>(9) * 10 + 111);
  ;
  if (comm->getRank() > 1 && comm->getRank() < comm->getSize() - 2)
    TEST_EQUALITY(xmap->getGlobalElement(11), (comm->getRank() + 1) * 40 + Teuchos::as<GO>(0) * 10 + 111);

  TEST_EQUALITY(xmap->getMinAllGlobalIndex(), novlxMap->getMinAllGlobalIndex());
  TEST_EQUALITY(xmap->getMaxAllGlobalIndex(), novlxMap->getMaxAllGlobalIndex());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, Merge, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 5;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const MatrixClass> A2 = bop->Merge();

  Teuchos::RCP<const BlockedCrsMatrixClass> bop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A2);
  TEST_EQUALITY(bop2, Teuchos::null);

  Teuchos::RCP<VectorClass> v1 = VectorFactoryClass::Build(bop->getRangeMap(), true);
  Teuchos::RCP<VectorClass> v2 = VectorFactoryClass::Build(A2->getRangeMap(), true);
  bop->getLocalDiagCopy(*v1);
  A2->getLocalDiagCopy(*v2);

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  magnitudeType tol   = Teuchos::ScalarTraits<magnitudeType>::eps();

  v1->update(-Teuchos::ScalarTraits<Scalar>::one(), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEUCHOS_TEST_COMPARE(v1->norm2(), <, tol, out, success);
  TEUCHOS_TEST_COMPARE(v1->normInf(), <, tol, out, success);

  TEST_EQUALITY(bop->getLocalNumEntries(), A2->getLocalNumEntries());
  TEST_EQUALITY(bop->getGlobalNumEntries(), A2->getGlobalNumEntries());
  TEST_EQUALITY(bop->getFrobeniusNorm(), A2->getFrobeniusNorm());

  // Thyra blocked operator
  Teuchos::RCP<BlockedCrsMatrixClass> bop3 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop3->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop3->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const MatrixClass> A4             = bop3->Merge();
  Teuchos::RCP<const BlockedCrsMatrixClass> bop4 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A4);
  TEST_EQUALITY(bop4, Teuchos::null);

  v1 = VectorFactoryClass::Build(bop3->getRangeMap(), true);
  v2 = VectorFactoryClass::Build(A4->getRangeMap(), true);
  bop3->getLocalDiagCopy(*v1);
  A4->getLocalDiagCopy(*v2);

  v1->update(-Teuchos::ScalarTraits<Scalar>::one(), *v2, Teuchos::ScalarTraits<Scalar>::one());
  TEUCHOS_TEST_COMPARE(v1->norm2(), <, tol, out, success);
  TEUCHOS_TEST_COMPARE(v1->normInf(), <, tol, out, success);

  TEST_EQUALITY(bop3->getLocalNumEntries(), A4->getLocalNumEntries());
  TEST_EQUALITY(bop3->getGlobalNumEntries(), A4->getGlobalNumEntries());
  TEST_EQUALITY(bop3->getFrobeniusNorm(), A4->getFrobeniusNorm());

  // Nested Xpetra blocked operator
  bop                                                 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 1 [ [ 2 4 0 ] 3] ]");

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);

  Teuchos::RCP<const MatrixClass> A               = brop->Merge();
  Teuchos::RCP<const BlockedCrsMatrixClass> brop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A);
  TEST_EQUALITY(brop2, Teuchos::null);

  v1 = VectorFactoryClass::Build(brop->getRangeMap(), true);
  v2 = VectorFactoryClass::Build(A->getRangeMap(), true);
  brop->getLocalDiagCopy(*v1);
  A->getLocalDiagCopy(*v2);

  // note that v1 and v2 have a different map here!
  // TEST_EQUALITY(v1->norm2(), v2->norm2());
  TEST_EQUALITY(v1->normInf(), v2->normInf());
  TEST_EQUALITY(v1->getMap()->isSameAs(*(v2->getMap())), false);
  TEST_EQUALITY(brop->getFullRangeMap()->isSameAs(*(A->getRangeMap())), true);

  TEST_EQUALITY(bop->getLocalNumEntries(), A->getLocalNumEntries());
  TEST_EQUALITY(bop->getGlobalNumEntries(), A->getGlobalNumEntries());
  TEST_EQUALITY(bop->getFrobeniusNorm(), A->getFrobeniusNorm());

  // Nested Thyra blocked operator
  bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  brop = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);

  A     = brop->Merge();
  brop2 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(A);
  TEST_EQUALITY(brop2, Teuchos::null);

  v1 = VectorFactoryClass::Build(brop->getRangeMap(), true);
  v2 = VectorFactoryClass::Build(A->getRangeMap(), true);
  brop->getLocalDiagCopy(*v1);
  A->getLocalDiagCopy(*v2);

  // note that v1 and v2 have the same map in thyra mode!
  // TEST_EQUALITY(v1->norm2(), v2->norm2());
  TEST_EQUALITY(v1->normInf(), v2->normInf());
  TEST_EQUALITY(v1->getMap()->isSameAs(*(v2->getMap())), false);
  TEST_EQUALITY(brop->getFullRangeMap()->isSameAs(*(A->getRangeMap())), true);

  TEST_EQUALITY(brop->getLocalNumEntries(), A->getLocalNumEntries());
  TEST_EQUALITY(brop->getGlobalNumEntries(), A->getGlobalNumEntries());
  TEUCHOS_TEST_COMPARE(Teuchos::ScalarTraits<Scalar>::magnitude(brop->getFrobeniusNorm() - A->getFrobeniusNorm()), <, 1e3 * tol, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, MatrixMatrixAdd, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 3;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const BlockedCrsMatrixClass> bop2 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop2->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop2->Cols(), Teuchos::as<size_t>(noBlocks));

  // matrix-matrix multiplication of blocked operators
  // Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar,LO,GO,Node> > bOpbOp = Xpetra::MatrixMatrix<Scalar,LO,GO,Node>::TwoMatrixMultiplyBlock(*bop,false,*bop2,false,out);
  Teuchos::RCP<MatrixClass> bOpOp = Teuchos::null;

  Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixAdd(*bop, false, STS::one(), *bop2, false, STS::one() * Teuchos::as<Scalar>(3.0), bOpOp, out);

  TEST_EQUALITY(bOpOp != Teuchos::null, true);
  Teuchos::RCP<BlockedCrsMatrixClass> bOpOp2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(bOpOp);
  TEST_EQUALITY(bOpOp2 != Teuchos::null, true);

  TEST_EQUALITY(bOpOp2->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bOpOp2->Cols(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor(), bop->getRangeMapExtractor());
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor(), bop2->getDomainMapExtractor());
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getFullMap()->isSameAs(*(bop->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getFullMap()->isSameAs(*(bop->getDomainMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getRangeMap()->isSameAs(*(bop->getMatrix(0, 1)->getRangeMap())), true);
  // TEST_EQUALITY(bOpOp2->getMatrix(0,1)->getDomainMap()->isSameAs(*(bop->getMatrix(0,1)->getDomainMap())),true);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getRangeMap()->isSameAs(*(bop->getMatrix(1, 0)->getRangeMap())), true);
  // TEST_EQUALITY(bOpOp2->getMatrix(1,0)->getDomainMap()->isSameAs(*(bop->getMatrix(1,0)->getDomainMap())),true);

  bOpOp2->fillComplete();

  TEST_EQUALITY(bOpOp2->getFrobeniusNorm(), 4.0 * bop->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getFrobeniusNorm(), STS::magnitude(STS::zero()));
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getFrobeniusNorm(), STS::magnitude(STS::zero()));
  TEST_EQUALITY(bOpOp2->getMatrix(1, 1)->getFrobeniusNorm(), 4.0 * bop->getMatrix(1, 1)->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(0, 0)->getFrobeniusNorm(), 4.0 * bop->getMatrix(0, 0)->getFrobeniusNorm());

  // Nested addition test (Xpetra)
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [1 2] ]");
  bop                                                 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  bop2                                                = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop2 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop2));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), false);

  bOpOp = Teuchos::null;

  Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixAdd(*brop, false, STS::one(), *brop2, false, STS::one() * Teuchos::as<Scalar>(3.0), bOpOp, out);

  TEST_EQUALITY(bOpOp != Teuchos::null, true);
  bOpOp2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(bOpOp);
  TEST_EQUALITY(bOpOp2 != Teuchos::null, true);

  TEST_EQUALITY(bOpOp2->Rows(), 2);
  TEST_EQUALITY(bOpOp2->Cols(), 2);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor(), brop->getRangeMapExtractor());
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor(), brop2->getDomainMapExtractor());
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getFullMap()->isSameAs(*(brop->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getFullMap()->isSameAs(*(brop->getDomainMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getRangeMap()->isSameAs(*(brop->getMatrix(0, 1)->getRangeMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getDomainMap()->isSameAs(*(brop->getMatrix(0, 1)->getDomainMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getRangeMap()->isSameAs(*(brop->getMatrix(1, 0)->getRangeMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getDomainMap()->isSameAs(*(brop->getMatrix(1, 0)->getDomainMap())), true);

  bOpOp2->fillComplete();

  TEST_EQUALITY(bOpOp2->getFrobeniusNorm(), 4.0 * brop->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getFrobeniusNorm(), STS::magnitude(STS::zero()));
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getFrobeniusNorm(), STS::magnitude(STS::zero()));
  TEST_EQUALITY(bOpOp2->getMatrix(1, 1)->getFrobeniusNorm(), 4.0 * brop2->getMatrix(1, 1)->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(0, 0)->getFrobeniusNorm(), 4.0 * brop2->getMatrix(0, 0)->getFrobeniusNorm());

  // Nested addition test (Thyra)
  bop  = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  bop2 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  brop  = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  brop2 = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop2));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), true);

  bOpOp = Teuchos::null;

  Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixAdd(*brop, false, STS::one(), *brop2, false, STS::one() * Teuchos::as<Scalar>(3.0), bOpOp, out);

  TEST_EQUALITY(bOpOp != Teuchos::null, true);
  bOpOp2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(bOpOp);
  TEST_EQUALITY(bOpOp2 != Teuchos::null, true);

  TEST_EQUALITY(bOpOp2->Rows(), 2);
  TEST_EQUALITY(bOpOp2->Cols(), 2);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor(), brop->getRangeMapExtractor());
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor(), brop2->getDomainMapExtractor());
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getFullMap()->isSameAs(*(brop->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getFullMap()->isSameAs(*(brop->getDomainMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getRangeMap()->isSameAs(*(brop->getMatrix(0, 1)->getRangeMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getDomainMap()->isSameAs(*(brop->getMatrix(0, 1)->getDomainMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getRangeMap()->isSameAs(*(brop->getMatrix(1, 0)->getRangeMap())), true);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getDomainMap()->isSameAs(*(brop->getMatrix(1, 0)->getDomainMap())), true);

  bOpOp2->fillComplete();

  TEST_EQUALITY(bOpOp2->getFrobeniusNorm(), 4.0 * brop->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1)->getFrobeniusNorm(), STS::magnitude(STS::zero()));
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0)->getFrobeniusNorm(), STS::magnitude(STS::zero()));
  TEST_EQUALITY(bOpOp2->getMatrix(1, 1)->getFrobeniusNorm(), 4.0 * brop2->getMatrix(1, 1)->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(0, 0)->getFrobeniusNorm(), 4.0 * brop2->getMatrix(0, 0)->getFrobeniusNorm());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, MatrixMatrixMultDiag, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrixClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  // typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                                  = 3;
  Teuchos::RCP<const BlockedCrsMatrixClass> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop->Cols(), Teuchos::as<size_t>(noBlocks));

  Teuchos::RCP<const BlockedCrsMatrixClass> bop2 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  TEST_EQUALITY(bop2->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bop2->Cols(), Teuchos::as<size_t>(noBlocks));

  // matrix-matrix multiplication of blocked operators
  Teuchos::RCP<MatrixClass> bOpOp = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixMultiplyBlock(*bop, false, *bop2, false, out);

  TEST_EQUALITY(bOpOp != Teuchos::null, true);
  Teuchos::RCP<BlockedCrsMatrixClass> bOpOp2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(bOpOp);
  TEST_EQUALITY(bOpOp2 != Teuchos::null, true);

  TEST_EQUALITY(bOpOp2->Rows(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bOpOp2->Cols(), Teuchos::as<size_t>(noBlocks));
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor(), bop->getRangeMapExtractor());
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor(), bop2->getDomainMapExtractor());
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getFullMap()->isSameAs(*(bop->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getFullMap()->isSameAs(*(bop2->getDomainMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getMap(1)->isSameAs(*(bop->getRangeMapExtractor()->getMap(1))), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getMap(1)->isSameAs(*(bop->getDomainMapExtractor()->getMap(1))), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 2), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(2, 0), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(2, 1), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 2), Teuchos::null);

  bOpOp2->fillComplete();

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  magnitudeType tol   = 1e6 * Teuchos::ScalarTraits<magnitudeType>::eps();

  TEST_EQUALITY(bOpOp2->getMatrix(0, 0)->getFrobeniusNorm(), bop->getMatrix(0, 0)->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp2->getMatrix(1, 1)->getFrobeniusNorm(), 2.0 * bop->getMatrix(1, 1)->getFrobeniusNorm());
  TEST_COMPARE(bOpOp2->getMatrix(2, 2)->getFrobeniusNorm() - 3.0 * bop->getMatrix(2, 2)->getFrobeniusNorm(), <, tol);

  Teuchos::RCP<VectorClass> v = VectorFactoryClass::Build(bOpOp2->getMatrix(2, 2)->getRangeMap(), true);
  bOpOp2->getMatrix(2, 2)->getLocalDiagCopy(*v);

  Teuchos::ArrayRCP<const Scalar> vdata = v->getData(0);
  bool bCheck                           = true;
  for (int i = 0; i < 10; i++)
    if (vdata[i] != Teuchos::as<Scalar>(9.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // Nested addition test (Xpetra)
  bop  = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  bop2 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ 0 [1 2] ]");
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  Teuchos::RCP<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>> brop2 =
      Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop2));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), false);

  // matrix-matrix multiplication of blocked operators
  bOpOp = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixMultiplyBlock(*brop, false, *brop2, false, out);

  TEST_EQUALITY(bOpOp != Teuchos::null, true);
  bOpOp2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(bOpOp);
  TEST_EQUALITY(bOpOp2 != Teuchos::null, true);

  TEST_EQUALITY(bOpOp2->Rows(), 2);
  TEST_EQUALITY(bOpOp2->Cols(), 2);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor(), brop->getRangeMapExtractor());
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor(), brop2->getDomainMapExtractor());
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getFullMap()->isSameAs(*(brop->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getFullMap()->isSameAs(*(brop2->getDomainMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getMap(1)->isSameAs(*(brop->getRangeMapExtractor()->getMap(1))), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getMap(1)->isSameAs(*(brop->getDomainMapExtractor()->getMap(1))), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0), Teuchos::null);

  Teuchos::RCP<const BlockedCrsMatrixClass> bOpOp21 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(bOpOp2->getMatrix(1, 1));
  TEST_EQUALITY(bOpOp21->Rows(), 2);
  TEST_EQUALITY(bOpOp21->Cols(), 2);
  TEST_EQUALITY(bOpOp21->getMatrix(0, 1), Teuchos::null);
  TEST_EQUALITY(bOpOp21->getMatrix(1, 0), Teuchos::null);

  TEST_EQUALITY(bOpOp2->getMatrix(0, 0)->getFrobeniusNorm(), bop->getMatrix(0, 0)->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp21->getMatrix(0, 0)->getFrobeniusNorm(), 2.0 * bop->getMatrix(1, 1)->getFrobeniusNorm());
  TEST_COMPARE(bOpOp21->getMatrix(1, 1)->getFrobeniusNorm() - 3.0 * bop->getMatrix(2, 2)->getFrobeniusNorm(), <, tol);

  v = VectorFactoryClass::Build(bOpOp21->getMatrix(1, 1)->getRangeMap(), true);
  bOpOp21->getMatrix(1, 1)->getLocalDiagCopy(*v);

  vdata  = v->getData(0);
  bCheck = true;
  for (int i = 0; i < 10; i++)
    if (vdata[i] != Teuchos::as<Scalar>(9.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);

  // Nested addition test (Thyra)
  bop  = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);
  bop2 = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrixThyra<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  brop  = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop));
  brop2 = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LO, GO, Node>>(buildReorderedBlockedCrsMatrix(brm, bop2));

  TEST_EQUALITY(brop->Rows(), 2);
  TEST_EQUALITY(brop->Cols(), 2);
  TEST_EQUALITY(brop->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop->getDomainMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop2->Rows(), 2);
  TEST_EQUALITY(brop2->Cols(), 2);
  TEST_EQUALITY(brop2->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(brop2->getDomainMapExtractor()->getThyraMode(), true);

  // matrix-matrix multiplication of blocked operators
  bOpOp = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixMultiplyBlock(*brop, false, *brop2, false, out);

  TEST_EQUALITY(bOpOp != Teuchos::null, true);
  bOpOp2 = Teuchos::rcp_dynamic_cast<BlockedCrsMatrixClass>(bOpOp);
  TEST_EQUALITY(bOpOp2 != Teuchos::null, true);

  TEST_EQUALITY(bOpOp2->Rows(), 2);
  TEST_EQUALITY(bOpOp2->Cols(), 2);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor(), brop->getRangeMapExtractor());
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor(), brop2->getDomainMapExtractor());
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getFullMap()->isSameAs(*(brop->getRangeMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getFullMap()->isSameAs(*(brop2->getDomainMapExtractor()->getFullMap())), true);
  TEST_EQUALITY(bOpOp2->getRangeMapExtractor()->getMap(1)->isSameAs(*(brop->getRangeMapExtractor()->getMap(1))), true);
  TEST_EQUALITY(bOpOp2->getDomainMapExtractor()->getMap(1)->isSameAs(*(brop->getDomainMapExtractor()->getMap(1))), true);
  TEST_EQUALITY(bOpOp2->getMatrix(0, 1), Teuchos::null);
  TEST_EQUALITY(bOpOp2->getMatrix(1, 0), Teuchos::null);

  bOpOp21 = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrixClass>(bOpOp2->getMatrix(1, 1));
  TEST_EQUALITY(bOpOp21->Rows(), 2);
  TEST_EQUALITY(bOpOp21->Cols(), 2);
  TEST_EQUALITY(bOpOp21->getMatrix(0, 1), Teuchos::null);
  TEST_EQUALITY(bOpOp21->getMatrix(1, 0), Teuchos::null);

  TEST_EQUALITY(bOpOp2->getMatrix(0, 0)->getFrobeniusNorm(), bop->getMatrix(0, 0)->getFrobeniusNorm());
  TEST_EQUALITY(bOpOp21->getMatrix(0, 0)->getFrobeniusNorm(), 2.0 * bop->getMatrix(1, 1)->getFrobeniusNorm());
  TEST_COMPARE(bOpOp21->getMatrix(1, 1)->getFrobeniusNorm() - 3.0 * bop->getMatrix(2, 2)->getFrobeniusNorm(), <, tol);

  v = VectorFactoryClass::Build(bOpOp21->getMatrix(1, 1)->getRangeMap(), true);
  bOpOp21->getMatrix(1, 1)->getLocalDiagCopy(*v);

  vdata  = v->getData(0);
  bCheck = true;
  for (int i = 0; i < 10; i++)
    if (vdata[i] != Teuchos::as<Scalar>(9.0)) bCheck = false;
  TEST_EQUALITY(bCheck, true);
}

/// simple test routine for the apply function of BlockedCrsMatrix
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, MatrixMatrixMult, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> MapClass;
  typedef Xpetra::StridedMap<LO, GO, Node> StridedMapClass;
  typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;
  typedef Xpetra::StridedMapFactory<LO, GO, Node> StridedMapFactoryClass;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> VectorClass;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> VectorFactoryClass;
  typedef Xpetra::MapExtractor<Scalar, LO, GO, Node> MapExtractorClass;
  typedef Xpetra::MapExtractorFactory<Scalar, LO, GO, Node> MapExtractorFactoryClass;
  typedef Xpetra::Matrix<Scalar, LO, GO, Node> MatrixClass;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const Teuchos::RCP<const MapClass> pointmap = MapFactoryClass::Build(lib, 12, 0, comm);

  // generate local maps for loading matrices
  Teuchos::Array<GO> velgidvec;  // global strided maps
  Teuchos::Array<GO> pregidvec;
  Teuchos::Array<GO> fullgidvec;  // full global map
  for (LO i = 0; i < Teuchos::as<LO>(pointmap->getLocalNumElements()); i++) {
    // loop over all local ids in pointmap

    // get corresponding global id
    GO gid = pointmap->getGlobalElement(i);

    // store global strided gids
    velgidvec.push_back(3 * gid);
    velgidvec.push_back(3 * gid + 1);
    pregidvec.push_back(3 * gid + 2);

    // gid for full map
    fullgidvec.push_back(3 * gid);
    fullgidvec.push_back(3 * gid + 1);
    fullgidvec.push_back(3 * gid + 2);
  }

  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  const Teuchos::RCP<const StridedMapClass> velmap  = StridedMapFactoryClass::Build(lib,
                                                                                    Teuchos::OrdinalTraits<GO>::invalid(),
                                                                                    velgidvec(),
                                                                                    0,
                                                                                    stridingInfo,
                                                                                    comm, 0);
  const Teuchos::RCP<const StridedMapClass> premap  = StridedMapFactoryClass::Build(lib,
                                                                                    Teuchos::OrdinalTraits<GO>::invalid(),
                                                                                    pregidvec(),
                                                                                    0,
                                                                                    stridingInfo,
                                                                                    comm, 1);
  const Teuchos::RCP<const StridedMapClass> fullmap = StridedMapFactoryClass::Build(lib,
                                                                                    Teuchos::OrdinalTraits<GO>::invalid(),
                                                                                    fullgidvec(),
                                                                                    0,
                                                                                    stridingInfo,
                                                                                    comm, -1);

  std::string tname = typeid(Scalar).name();
  if (tname.find("complex") != std::string::npos) {
    std::cout << "Skip test for scalar=" << tname << std::endl;
    return;
  }

  Teuchos::RCP<MatrixClass> A = Xpetra::IO<Scalar, LO, GO, Node>::Read("A.mat", fullmap->getMap());

  std::vector<Teuchos::RCP<const MapClass>> xmaps;
  xmaps.push_back(velmap);
  xmaps.push_back(premap);

  Teuchos::RCP<const MapExtractorClass> rowMapExtractormap_extractor = MapExtractorFactoryClass::Build(fullmap, xmaps);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bOp =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*A, rowMapExtractormap_extractor, rowMapExtractormap_extractor);

  Teuchos::RCP<MatrixClass> A2 = Xpetra::IO<Scalar, LO, GO, Node>::Read("A.mat", fullmap->getMap());

  std::vector<Teuchos::RCP<const MapClass>> xmaps2;
  xmaps2.push_back(velmap);
  xmaps2.push_back(premap);

  Teuchos::RCP<const MapExtractorClass> map_extractor2 = MapExtractorFactoryClass::Build(fullmap, xmaps2);

  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bOp2 =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*A2, map_extractor2, map_extractor2);

  // matrix-matrix multiplication of standard matrices
  Teuchos::RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> fuAfuA_2 = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply(*A, false, *A2, false, out);
  fuAfuA_2->describe(out);

  // matrix-matrix multiplication of blocked operators
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> bOpbOp_2 = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixMultiplyBlock(*bOp, false, *bOp2, false, out);

  // build gloabl vector with one entries
  Teuchos::RCP<VectorClass> ones = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> exp  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> res  = VectorFactoryClass::Build(fullmap, true);
  Teuchos::RCP<VectorClass> rnd  = VectorFactoryClass::Build(fullmap, true);
  ones->putScalar(STS::one());
  rnd->randomize();

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  magnitudeType tol1  = Teuchos::ScalarTraits<magnitudeType>::eps();
  magnitudeType tol2  = 500 * tol1;

  fuAfuA_2->apply(*ones, *exp);
  bOpbOp_2->apply(*ones, *res);
  res->update(-STS::one(), *exp, STS::one());
  TEUCHOS_TEST_COMPARE(res->norm2(), <, tol1, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, tol1, out, success);

  A->apply(*rnd, *exp);
  bOp->apply(*rnd, *res);
  res->update(-STS::one(), *exp, STS::one());

  TEUCHOS_TEST_COMPARE(res->norm2(), <, tol2, out, success);
  TEUCHOS_TEST_COMPARE(res->normInf(), <, tol2, out, success);

  TEUCHOS_TEST_EQUALITY(fuAfuA_2->getGlobalNumEntries(), 312, out, success);
  TEUCHOS_TEST_EQUALITY(bOpbOp_2->getGlobalNumEntries(), 312, out, success);

  Teuchos::RCP<const MapClass> rgMap0           = bOpbOp_2->getRangeMap(0);
  Teuchos::RCP<const StridedMapClass> strRgMap0 = Teuchos::rcp_dynamic_cast<const StridedMapClass>(rgMap0);
  TEUCHOS_TEST_EQUALITY(strRgMap0 == Teuchos::null, false, out, success);
  std::vector<size_t> strInfoData = strRgMap0->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success);
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success);
  TEUCHOS_TEST_EQUALITY(strRgMap0->getFixedBlockSize(), 3, out, success);
  TEUCHOS_TEST_EQUALITY(strRgMap0->getStridedBlockId(), 0, out, success);

  /* TODO think about this
  Teuchos::RCP<const MapClass> rgMap = bOpbOp_2->getRangeMap();
  Teuchos::RCP<const StridedMapClass> strRgMap = Teuchos::rcp_dynamic_cast<const StridedMapClass>(rgMap);
  TEUCHOS_TEST_EQUALITY(strRgMap==Teuchos::null, false, out, success );
  strInfoData = strRgMap->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success );
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success );
  TEUCHOS_TEST_EQUALITY(strRgMap->getFixedBlockSize(), 3, out, success );
  TEUCHOS_TEST_EQUALITY(strRgMap->getStridedBlockId(), -1, out, success );*/

  Teuchos::RCP<const MapClass> doMap0           = bOpbOp_2->getDomainMap(0);
  Teuchos::RCP<const StridedMapClass> strDoMap0 = Teuchos::rcp_dynamic_cast<const StridedMapClass>(doMap0);
  TEUCHOS_TEST_EQUALITY(strDoMap0 == Teuchos::null, false, out, success);
  strInfoData = strDoMap0->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success);
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success);
  TEUCHOS_TEST_EQUALITY(strDoMap0->getFixedBlockSize(), 3, out, success);
  TEUCHOS_TEST_EQUALITY(strDoMap0->getStridedBlockId(), 0, out, success);

  /* TODO think about this
  Teuchos::RCP<const MapClass> doMap = bOpbOp_2->getDomainMap();
  Teuchos::RCP<const StridedMapClass> strDoMap = Teuchos::rcp_dynamic_cast<const StridedMapClass>(doMap);
  TEUCHOS_TEST_EQUALITY(strDoMap==Teuchos::null, false, out, success );
  strInfoData = strDoMap->getStridingData();
  TEUCHOS_TEST_EQUALITY(strInfoData[0], 2, out, success );
  TEUCHOS_TEST_EQUALITY(strInfoData[1], 1, out, success );
  TEUCHOS_TEST_EQUALITY(strDoMap->getFixedBlockSize(), 3, out, success );
  TEUCHOS_TEST_EQUALITY(strDoMap->getStridedBlockId(), -1, out, success );
  */
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, BlockedOperatorApply, M, MA, Scalar, LO, GO, Node) {
  typedef Xpetra::Map<LO, GO, Node> Map;
  typedef Xpetra::BlockedMap<LO, GO, Node> BlockedMap;
  typedef Xpetra::MapExtractor<Scalar, LO, GO, Node> MapExtractor;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LO, GO, Node> BlockedMultiVector;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> MultiVectorFactory;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> BlockedCrsMatrix;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // get a comm and node
  Teuchos::RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  int noBlocks                             = 3;
  Teuchos::RCP<const BlockedCrsMatrix> bop = XpetraBlockMatrixTests::CreateBlockDiagonalExampleMatrix<Scalar, LO, GO, Node, M>(noBlocks, *comm);

  TEST_EQUALITY(bop->Rows(), 3);
  TEST_EQUALITY(bop->Cols(), 3);
  TEST_EQUALITY(bop->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bop->getDomainMapExtractor()->getThyraMode(), false);

  // build gloabl vector with one entries (build monolithic maps)
  Teuchos::RCP<const Map> rgMap         = bop->getRangeMap();
  Teuchos::RCP<const BlockedMap> rgBMap = Teuchos::rcp_dynamic_cast<const BlockedMap>(rgMap);
  Teuchos::RCP<MultiVector> ones        = MultiVectorFactory::Build(rgBMap->getFullMap(), 1, true);
  Teuchos::RCP<MultiVector> res1        = MultiVectorFactory::Build(rgBMap->getFullMap(), 1, true);
  ones->putScalar(STS::one());
  res1->putScalar(STS::zero());

  Teuchos::RCP<const MapExtractor> meRange  = bop->getRangeMapExtractor();
  Teuchos::RCP<const MapExtractor> meDomain = bop->getDomainMapExtractor();

  // create BlockedMultiVectors
  Teuchos::RCP<BlockedMultiVector> bones =
      Teuchos::rcp(new BlockedMultiVector(meDomain, ones));
  Teuchos::RCP<BlockedMultiVector> res2 =
      Teuchos::rcp(new BlockedMultiVector(meRange, res1));
  Teuchos::RCP<BlockedMultiVector> res3 =
      Teuchos::rcp(new BlockedMultiVector(meRange, res1));

  // input blocked, output standard
  TEST_NOTHROW(bop->apply(*bones, *res1));
  // input blocked, output blocked
  TEST_NOTHROW(bop->apply(*bones, *res2));
  // input standard, output blocked
  TEST_NOTHROW(bop->apply(*ones, *res3));

  for (size_t r = 0; r < meRange->NumMaps(); r++) {
    Teuchos::RCP<MultiVector> part2 = meRange->ExtractVector(res2, r);
    Teuchos::RCP<MultiVector> part3 = meRange->ExtractVector(res3, r);

    Teuchos::ArrayRCP<const Scalar> partd2 = part2->getData(0);
    Teuchos::ArrayRCP<const Scalar> partd3 = part3->getData(0);
    for (LO l = 0; l < Teuchos::as<LO>(part2->getLocalLength()); l++) {
      TEST_EQUALITY(partd2[l], Teuchos::as<Scalar>(r + 1) * STS::one());
      TEST_EQUALITY(partd3[l], Teuchos::as<Scalar>(r + 1) * STS::one());
    }
  }

  Teuchos::RCP<MultiVector> merged_res2 = res2->Merge();
  Teuchos::ArrayRCP<const Scalar> resd1 = res1->getData(0);
  Teuchos::ArrayRCP<const Scalar> resd2 = merged_res2->getData(0);

  for (LO l = 0; l < Teuchos::as<LO>(res1->getLocalLength()); l++) {
    TEST_EQUALITY(resd1[l], resd2[l]);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ConstructFromBlockedMap, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using BlockedCrsMatrix = Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>;
  using BlockedMap       = Xpetra::BlockedMap<LO, GO, Node>;
  using Map              = Xpetra::Map<LO, GO, Node>;
  using MapFactory       = Xpetra::MapFactory<LO, GO, Node>;
  using MapUtils         = Xpetra::MapUtils<LO, GO, Node>;

  RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const GO gNumElementsPerBlock = 29;
  Array<GO> gidsOne;

  RCP<const Map> map1 = MapFactory::Build(lib, gNumElementsPerBlock, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_ASSERT(!map1.is_null());
  TEST_EQUALITY_CONST(map1->getGlobalNumElements(), gNumElementsPerBlock);

  ArrayView<const GO> myGIDs1 = map1->getLocalElementList();
  Array<GO> myGIDs2;
  for (const auto& gid1 : myGIDs1)
    myGIDs2.push_back(gid1 + gNumElementsPerBlock);
  RCP<const Map> map2 = MapFactory::Build(lib, gNumElementsPerBlock, myGIDs2, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_ASSERT(!map2.is_null());
  TEST_EQUALITY_CONST(map2->getGlobalNumElements(), gNumElementsPerBlock);

  std::vector<RCP<const Map>> maps;
  maps.push_back(map1);
  maps.push_back(map2);
  RCP<const Map> fullMap           = MapUtils::concatenateMaps(maps);
  RCP<const BlockedMap> blockedMap = rcp(new BlockedMap(fullMap, maps));

  TEST_ASSERT(!blockedMap.is_null());
  TEST_EQUALITY(blockedMap->getNumMaps(), 2);

  RCP<BlockedCrsMatrix> blockMatrix = rcp(new BlockedCrsMatrix(blockedMap, blockedMap, 1));
  TEST_ASSERT(!blockMatrix.is_null());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ConstructFromMapExtractor, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using BlockedCrsMatrix    = Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>;
  using Map                 = Xpetra::Map<LO, GO, Node>;
  using MapExtractor        = Xpetra::MapExtractor<Scalar, LO, GO, Node>;
  using MapExtractorFactory = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>;
  using MapFactory          = Xpetra::MapFactory<LO, GO, Node>;
  using MapUtils            = Xpetra::MapUtils<LO, GO, Node>;

  RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const GO gNumElementsPerBlock = 29;
  Array<GO> gidsOne;

  RCP<const Map> map1 = MapFactory::Build(lib, gNumElementsPerBlock, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_ASSERT(!map1.is_null());
  TEST_EQUALITY_CONST(map1->getGlobalNumElements(), gNumElementsPerBlock);

  ArrayView<const GO> myGIDs1 = map1->getLocalElementList();
  Array<GO> myGIDs2;
  for (const auto& gid1 : myGIDs1)
    myGIDs2.push_back(gid1 + gNumElementsPerBlock);
  RCP<const Map> map2 = MapFactory::Build(lib, gNumElementsPerBlock, myGIDs2, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_ASSERT(!map2.is_null());
  TEST_EQUALITY_CONST(map2->getGlobalNumElements(), gNumElementsPerBlock);

  std::vector<RCP<const Map>> maps;
  maps.push_back(map1);
  maps.push_back(map2);
  RCP<const Map> fullMap               = MapUtils::concatenateMaps(maps);
  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(fullMap, maps);

  TEST_ASSERT(!mapExtractor.is_null());

  RCP<BlockedCrsMatrix> blockMatrix = rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 1));
  TEST_ASSERT(!blockMatrix.is_null());
}

TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, ConstructFromBlockedVector, M, MA, Scalar, LO, GO, Node) {
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using BlockedCrsMatrix = Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>;
  using BlockedMap       = Xpetra::BlockedMap<LO, GO, Node>;
  using Map              = Xpetra::Map<LO, GO, Node>;
  using MapFactory       = Xpetra::MapFactory<LO, GO, Node>;
  using MapUtils         = Xpetra::MapUtils<LO, GO, Node>;

  using STS = Teuchos::ScalarTraits<Scalar>;
  typedef typename STS::magnitudeType MT;

  RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  M testMap(1, 0, comm);
  Xpetra::UnderlyingLib lib = testMap.lib();

  const GO gNumElementsPerBlock = 29;
  Array<GO> gidsOne;

  RCP<const Map> map1 = MapFactory::Build(lib, gNumElementsPerBlock, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_ASSERT(!map1.is_null());
  TEST_EQUALITY_CONST(map1->getGlobalNumElements(), gNumElementsPerBlock);

  ArrayView<const GO> myGIDs1 = map1->getLocalElementList();
  Array<GO> myGIDs2;
  for (const auto& gid1 : myGIDs1)
    myGIDs2.push_back(gid1 + gNumElementsPerBlock);
  RCP<const Map> map2 = MapFactory::Build(lib, gNumElementsPerBlock, myGIDs2, Teuchos::ScalarTraits<GO>::zero(), comm);

  TEST_ASSERT(!map2.is_null());
  TEST_EQUALITY_CONST(map2->getGlobalNumElements(), gNumElementsPerBlock);

  std::vector<RCP<const Map>> maps;
  maps.push_back(map1);
  maps.push_back(map2);
  RCP<const Map> fullMap           = MapUtils::concatenateMaps(maps);
  RCP<const BlockedMap> blockedMap = rcp(new BlockedMap(fullMap, maps));

  TEST_ASSERT(!blockedMap.is_null());
  TEST_EQUALITY(blockedMap->getNumMaps(), 2);

  const RCP<Xpetra::Vector<Scalar, LO, GO, Node>> vec = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(blockedMap);
  vec->randomize();

  RCP<Xpetra::BlockedVector<Scalar, LO, GO, Node>> blockVec =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedVector<Scalar, LO, GO, Node>>(vec);
  TEST_ASSERT(!blockVec.is_null());

  RCP<Xpetra::Matrix<Scalar, LO, GO, Node>> matrix = Xpetra::MatrixFactory<Scalar, LO, GO, Node>::Build(vec.getConst());
  TEST_ASSERT(!matrix.is_null());

  RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>> blockMatrix =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>>(matrix);
  TEST_ASSERT(!blockMatrix.is_null());

  const RCP<Xpetra::Vector<Scalar, LO, GO, Node>> diagonal = Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(blockedMap);
  blockMatrix->getLocalDiagCopy(*diagonal);
  TEST_ASSERT(!diagonal.is_null());

  RCP<Xpetra::BlockedVector<Scalar, LO, GO, Node>> blockDiagonal =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedVector<Scalar, LO, GO, Node>>(diagonal);
  TEST_ASSERT(!blockDiagonal.is_null());

  const MT tol = 1e-12;

  TEST_EQUALITY(gNumElementsPerBlock, blockMatrix->getMatrix(0, 0)->getGlobalNumEntries());
  TEST_FLOATING_EQUALITY(blockVec->getMultiVector(0)->getVector(0)->norm2(), blockDiagonal->getMultiVector(0)->getVector(0)->norm2(), tol);
  TEST_FLOATING_EQUALITY(blockVec->getMultiVector(1)->getVector(0)->norm2(), blockDiagonal->getMultiVector(1)->getVector(0)->norm2(), tol);
  TEST_EQUALITY(gNumElementsPerBlock, blockMatrix->getMatrix(1, 1)->getGlobalNumEntries());
  TEST_FLOATING_EQUALITY(blockVec->getMultiVector(0)->getVector(0)->norm2(), blockMatrix->getMatrix(0, 0)->getFrobeniusNorm(), tol);
  TEST_FLOATING_EQUALITY(blockVec->getMultiVector(1)->getVector(0)->norm2(), blockMatrix->getMatrix(1, 1)->getFrobeniusNorm(), tol);
  TEST_FLOATING_EQUALITY(blockVec->norm2(), blockDiagonal->norm2(), tol);
  TEST_FLOATING_EQUALITY(blockVec->norm2(), blockMatrix->getFrobeniusNorm(), tol);
}

// simple test for matrix-matrix multiplication for a 2x2 blocked matrix with a 2x1 blocked matrix
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(BlockedCrsMatrix, EpetraMatrixMatrixMult2x1, M, MA, Scalar, LO, GO, Node) {
#ifdef HAVE_XPETRA_EPETRAEXT
  using Teuchos::RCP;
  using Teuchos::rcp;

  using BlockedCrsMatrix    = Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>;
  using CrsMatrix           = Xpetra::CrsMatrix<Scalar, LO, GO, Node>;
  using CrsMatrixWrap       = Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>;
  using MapExtractor        = Xpetra::MapExtractor<Scalar, LO, GO, Node>;
  using MapExtractorFactory = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>;
  using Matrix              = Xpetra::Matrix<Scalar, LO, GO, Node>;
  using Vector              = Xpetra::Vector<Scalar, LO, GO, Node>;
  using VectorFactory       = Xpetra::VectorFactory<Scalar, LO, GO, Node>;

  RCP<const Teuchos::Comm<int>> comm = getDefaultComm();

  // build maps
  RCP<Epetra_Map> rowMap1    = rcp(new Epetra_Map(24, 0, *Xpetra::toEpetra(comm)));
  RCP<Epetra_Map> rowMap2    = rcp(new Epetra_Map(12, 24, *Xpetra::toEpetra(comm)));
  RCP<Epetra_Map> domainMap1 = rcp(new Epetra_Map(8, 0, *Xpetra::toEpetra(comm)));
  RCP<Epetra_Map> domainMap2 = rcp(new Epetra_Map(4, 8, *Xpetra::toEpetra(comm)));

  std::vector<RCP<const Epetra_Map>> rowMaps;
  rowMaps.push_back(rowMap1);
  rowMaps.push_back(rowMap2);
  std::vector<RCP<const Epetra_Map>> domainMaps;
  domainMaps.push_back(domainMap1);
  domainMaps.push_back(domainMap2);

  RCP<Epetra_Map> fullRowMap    = MergeMaps(rowMaps);
  RCP<Epetra_Map> fullDomainMap = MergeMaps(domainMaps);

  // read in matrices in matrix market format
  Epetra_CrsMatrix* ptrA = 0;
  Epetra_CrsMatrix* ptrP = 0;
  EpetraExt::MatrixMarketFileToCrsMatrix("A.mat", *fullRowMap, *fullRowMap, *fullRowMap, ptrA);
  EpetraExt::MatrixMarketFileToCrsMatrix("P.mat", *fullRowMap, *fullRowMap, *fullDomainMap, ptrP);
  RCP<Epetra_CrsMatrix> epA = rcp(ptrA);
  RCP<Epetra_CrsMatrix> epP = rcp(ptrP);

  // Transform Epetra stuff to Xpetra

  RCP<Xpetra::EpetraMapT<GO, Node>> xFullRowMap    = rcp(new Xpetra::EpetraMapT<GO, Node>(fullRowMap));
  RCP<Xpetra::EpetraMapT<GO, Node>> xFullDomainMap = rcp(new Xpetra::EpetraMapT<GO, Node>(fullDomainMap));
  RCP<Xpetra::EpetraMapT<GO, Node>> xRowMap1       = rcp(new Xpetra::EpetraMapT<GO, Node>(rowMap1));
  RCP<Xpetra::EpetraMapT<GO, Node>> xRowMap2       = rcp(new Xpetra::EpetraMapT<GO, Node>(rowMap2));
  RCP<Xpetra::EpetraMapT<GO, Node>> xDomainMap1    = rcp(new Xpetra::EpetraMapT<GO, Node>(domainMap1));
  RCP<Xpetra::EpetraMapT<GO, Node>> xDomainMap2    = rcp(new Xpetra::EpetraMapT<GO, Node>(domainMap2));

  // build map extractor objects
  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO>>> xRowMaps;
  xRowMaps.push_back(xRowMap1);
  xRowMaps.push_back(xRowMap2);
  RCP<const MapExtractor> rowMapExtractor = MapExtractorFactory::Build(xFullRowMap, xRowMaps);

  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO>>> xDomainMaps;
  xDomainMaps.push_back(xFullDomainMap);
  RCP<const MapExtractor> domainMapExtractor = MapExtractorFactory::Build(xFullDomainMap, xDomainMaps);

  // build blocked operators

  // build 2x2 blocked operator
  RCP<CrsMatrix> xCrsA = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(epA));
  RCP<Matrix> xA       = rcp(new CrsMatrixWrap(xCrsA));
  Teuchos::RCP<BlockedCrsMatrix> bA =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*xA, rowMapExtractor, rowMapExtractor);

  TEUCHOS_TEST_EQUALITY(bA->Rows(), 2, out, success);
  TEUCHOS_TEST_EQUALITY(bA->Cols(), 2, out, success);

  // build 2x1 blocked operator
  RCP<CrsMatrix> xCrsP = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(epP));
  RCP<Matrix> xP       = rcp(new CrsMatrixWrap(xCrsP));
  Teuchos::RCP<BlockedCrsMatrix> bP =
      Xpetra::MatrixUtils<Scalar, LO, GO, Node>::SplitMatrix(*xP, rowMapExtractor, domainMapExtractor);

  TEUCHOS_TEST_EQUALITY(bP->Rows(), 2, out, success);
  TEUCHOS_TEST_EQUALITY(bP->Cols(), 1, out, success);

  RCP<BlockedCrsMatrix> bAbP = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::TwoMatrixMultiplyBlock(*bA, false, *bP, false, out);

  TEUCHOS_TEST_EQUALITY(bAbP->Rows(), 2, out, success);
  TEUCHOS_TEST_EQUALITY(bAbP->Cols(), 1, out, success);

  RCP<Matrix> xAP = Xpetra::MatrixMatrix<Scalar, LO, GO, Node>::Multiply(*xA, false, *xP, false, out);

  // Test if blocked and merged MatVec deliver the same result
  RCP<Vector> oneVectorBlocked = VectorFactory::Build(bAbP->getDomainMap(), true);
  RCP<Vector> resVectorBlocked = VectorFactory::Build(bAbP->getRangeMap(), true);
  oneVectorBlocked->putScalar(Teuchos::ScalarTraits<Scalar>::one());
  bAbP->apply(*oneVectorBlocked, *resVectorBlocked);
  TEUCHOS_TEST_COMPARE(resVectorBlocked->norm2(), >, 1.0e-16, out, success);

  RCP<Vector> oneVector = VectorFactory::Build(xAP->getDomainMap(), true);
  RCP<Vector> resVector = VectorFactory::Build(xAP->getRangeMap(), true);
  oneVector->putScalar(Teuchos::ScalarTraits<Scalar>::one());
  xAP->apply(*oneVector, *resVector);
  TEUCHOS_TEST_COMPARE(resVector->norm2(), >, 1.0e-16, out, success);

  resVectorBlocked->update(1.0, *resVector, -1.0);
  TEUCHOS_TEST_COMPARE(resVectorBlocked->normInf(), <, 1.0e-16, out, success);
#endif
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

#define XP_MATRIX_INSTANT(S, LO, GO, N)                                                                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, SplitMatrix, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperator, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperator2, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperatorThyra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperator2Thyra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperatorApply, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperatorApply2, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperatorApplyThyra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperatorApplyThyraSmall, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReorderBlockOperatorApply2Thyra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, Apply, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, getLocalDiagCopy, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, leftScale, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, rightScale, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, deepCopy, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, transformThyra2XpetraGIDs, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, Merge, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, MatrixMatrixAdd, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, MatrixMatrixMultDiag, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, MatrixMatrixMult, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, BlockedOperatorApply, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ConstructFromBlockedMap, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ConstructFromMapExtractor, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ConstructFromBlockedVector, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

// List of tests which run only with Tpetra
#define XP_TPETRA_MATRIX_INSTANT(S, LO, GO, N)                                                                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, CreateBlockedDiagonalOp, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, CreateBlockedDiagonalOpThyra, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

// List of tests which run only with Epetra
#define XP_EPETRA_MATRIX_INSTANT(S, LO, GO, N)                                                                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, ReadWriteBlockedMatrix, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(BlockedCrsMatrix, EpetraMatrixMatrixMult2x1, M##LO##GO##N, MA##S##LO##GO##N, S, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_TPETRA_MATRIX_INSTANT)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_MATRIX_INSTANT)

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double, int, int, EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double, int, int, EpetraNode)
XP_MATRIX_INSTANT(double, int, int, EpetraNode)
#endif
// EpetraExt routines are not working with 64 bit
/*#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_TYPES(double,int,LongLong,EpetraNode)
XP_EPETRA_MATRIX_INSTANT(double,int,LongLong,EpetraNode)
#endif*/

#endif

}  // namespace XpetraBlockMatrixTests
