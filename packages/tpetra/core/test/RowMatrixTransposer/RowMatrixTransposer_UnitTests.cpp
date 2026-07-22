// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include "KokkosKernels_ArithTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ScalarTraits;

using Tpetra::createUniformContigMapWithNode;
using Tpetra::CrsMatrix;
using Tpetra::RowMatrixTransposer;
using Tpetra::MatrixMarket::Reader;
using Tpetra::MatrixMatrix::Add;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RowMatrixTransposer, RectangularTranspose, LO, GO, Node) {
  typedef CrsMatrix<>::scalar_type Scalar;
  typedef CrsMatrix<Scalar, LO, GO, Node> MAT;
  auto comm = Tpetra::getDefaultComm();

  RCP<MAT> matrix  = Reader<MAT>::readSparseFile("a.mtx", comm);
  RCP<MAT> matrixT = Reader<MAT>::readSparseFile("atrans.mtx", comm);

  RowMatrixTransposer<Scalar, LO, GO, Node> at(matrix);
  RCP<MAT> calculated = at.createTranspose();

  RCP<MAT> diffMatrix = rcp(new MAT(matrixT->getRowMap(), matrixT->getLocalMaxNumRowEntries()));

  Scalar sOne = ScalarTraits<Scalar>::one();
  Add(*calculated, false, -sOne, *matrixT, false, sOne, diffMatrix);
  diffMatrix->fillComplete(matrixT->getDomainMap(), matrixT->getRangeMap());

  Scalar diffNorm = diffMatrix->getFrobeniusNorm();
  Scalar realNorm = matrixT->getFrobeniusNorm();
  Scalar epsilon  = diffNorm / realNorm;

  TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon), <, 1e-10)
  TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon), <, 1e-10)
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(RowMatrixTransposer, ComplexRectangularTranspose, LO, GO, Node) {
  typedef CrsMatrix<>::scalar_type Scalar;
  typedef CrsMatrix<Scalar, LO, GO, Node> MAT;
#if KOKKOSKERNELS_VERSION < 50199
  // KokkosKernels does not yet support conjugating the values.
  return;
#endif
  if (!Teuchos::ScalarTraits<Scalar>::isComplex)
    return;
  auto comm = Tpetra::getDefaultComm();

  RCP<MAT> matrix  = Reader<MAT>::readSparseFile("a.mtx", comm);
  RCP<MAT> matrixT = Reader<MAT>::readSparseFile("atrans.mtx", comm);

  {
    // matrix is real valued. Let's make it imaginary.
    using KAT      = KokkosKernels::ArithTraits<typename MAT::impl_scalar_type>;
    auto lclMatrix = matrix->getLocalMatrixHost();
    auto values    = lclMatrix.values;
    auto j         = KAT::sqrt(-KAT::one());
    for (size_t i = 0; i < values.extent(0); ++i)
      values(i) = j * values(j);
  }
  {
    using KAT      = KokkosKernels::ArithTraits<typename MAT::impl_scalar_type>;
    auto lclMatrix = matrixT->getLocalMatrixHost();
    auto values    = lclMatrix.values;
    auto j         = KAT::sqrt(-KAT::one());
    for (size_t i = 0; i < values.extent(0); ++i)
      values(i) = KAT::conj(j * values(j));
  }

  RowMatrixTransposer<Scalar, LO, GO, Node> at(matrix);
  RCP<MAT> calculated = at.createTranspose();

  RCP<MAT> diffMatrix = rcp(new MAT(matrixT->getRowMap(), matrixT->getLocalMaxNumRowEntries()));

  Scalar sOne = ScalarTraits<Scalar>::one();
  Add(*calculated, false, -sOne, *matrixT, false, sOne, diffMatrix);
  diffMatrix->fillComplete(matrixT->getDomainMap(), matrixT->getRangeMap());

  Scalar diffNorm = diffMatrix->getFrobeniusNorm();
  Scalar realNorm = matrixT->getFrobeniusNorm();
  Scalar epsilon  = diffNorm / realNorm;

  TEST_COMPARE(ScalarTraits<Scalar>::real(epsilon), <, 1e-10)
  TEST_COMPARE(ScalarTraits<Scalar>::imag(epsilon), <, 1e-10)
}

#define UNIT_TEST_GROUP(LO, GO, NODE)                                                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(RowMatrixTransposer, RectangularTranspose, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(RowMatrixTransposer, ComplexRectangularTranspose, LO, GO, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

}  // namespace
