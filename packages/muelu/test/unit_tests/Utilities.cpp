// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Comm.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_BlockReorderManager.hpp>
#include <Xpetra_ReorderedBlockedCrsMatrix.hpp>
#include <Xpetra_IO.hpp>

#include <MueLu_config.hpp>
#include <MueLu_TestHelpers.hpp>
#include <MueLu_Version.hpp>
#include <MueLu_Utilities.hpp>

// This file is intended to house all the tests for MueLu_Utilities.hpp.

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, MatMatMult_EpetraVsTpetra, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  out << "version: " << MueLu::Version() << std::endl;
  out << "This test compares the matrix matrix multiply between Tpetra and Epetra" << std::endl;

  MUELU_TESTING_LIMIT_EPETRA_SCOPE_TPETRA_IS_DEFAULT(Scalar, GlobalOrdinal, Node);

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  // Calculate result = (Op*Op)*X for Epetra
  GO nx                   = 37 * comm->getSize();
  GO ny                   = nx;
  RCP<Matrix> Op          = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(nx, ny, Xpetra::UseEpetra);
  RCP<Matrix> OpOp        = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Op, false, *Op, false, out);
  RCP<MultiVector> result = MultiVectorFactory::Build(OpOp->getRangeMap(), 1);
  RCP<MultiVector> X      = MultiVectorFactory::Build(OpOp->getDomainMap(), 1);
  Teuchos::Array<magnitude_type> xnorm(1);
  X->setSeed(8675309);
  X->randomize(true);
  X->norm2(xnorm);
  OpOp->apply(*X, *result, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Teuchos::Array<magnitude_type> normEpetra(1);
  result->norm2(normEpetra);

  // aid debugging by calculating Op*(Op*X)
  RCP<MultiVector> workVec = MultiVectorFactory::Build(OpOp->getRangeMap(), 1);
  RCP<MultiVector> check1  = MultiVectorFactory::Build(OpOp->getRangeMap(), 1);
  Op->apply(*X, *workVec, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Op->apply(*workVec, *check1, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Teuchos::Array<magnitude_type> normCheck1(1);
  check1->norm2(normCheck1);

  // Calculate result = (Op*Op)*X for Tpetra
  Op     = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build2DPoisson(nx, ny, Xpetra::UseTpetra);
  OpOp   = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*Op, false, *Op, false, out);
  result = MultiVectorFactory::Build(OpOp->getRangeMap(), 1);
  X      = MultiVectorFactory::Build(OpOp->getDomainMap(), 1);
  X->setSeed(8675309);
  X->randomize(true);
  X->norm2(xnorm);
  OpOp->apply(*X, *result, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Teuchos::Array<magnitude_type> normTpetra(1);
  result->norm2(normTpetra);

  // aid debugging by calculating Op*(Op*X)
  workVec                 = MultiVectorFactory::Build(OpOp->getRangeMap(), 1);
  RCP<MultiVector> check2 = MultiVectorFactory::Build(OpOp->getRangeMap(), 1);
  Op->apply(*X, *workVec, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Op->apply(*workVec, *check2, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  Teuchos::Array<magnitude_type> normCheck2(1);
  check2->norm2(normCheck2);

  TEST_FLOATING_EQUALITY(normEpetra[0], normTpetra[0], 1e-12);
  out << "Epetra ||A*(A*x)|| = " << normCheck1[0] << std::endl;
  out << "Tpetra ||A*(A*x)|| = " << normCheck2[0] << std::endl;
#else
  out << "Skipping test because some required packages are not enabled (Tpetra, EpetraExt)." << std::endl;
#endif

}  // EpetraVersusTpetra

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, DetectDirichletRows, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build1DPoisson(100);
  Teuchos::ArrayView<const LocalOrdinal> indices;
  Teuchos::ArrayView<const Scalar> values;

  LocalOrdinal localRowToZero = 5;
  A->resumeFill();
  A->getLocalRowView(localRowToZero, indices, values);
  Array<Scalar> newvalues(values.size(), TST::zero());
  for (int j = 0; j < indices.size(); j++)
    // keep diagonal
    if (indices[j] == localRowToZero) newvalues[j] = values[j];
  A->replaceLocalValues(localRowToZero, indices, newvalues);

  A->fillComplete();

  ArrayRCP<const bool> drows = Utilities::DetectDirichletRows(*A);
  TEST_EQUALITY(drows[localRowToZero], true);
  TEST_EQUALITY(drows[localRowToZero - 1], false);

  A->resumeFill();
  A->getLocalRowView(localRowToZero, indices, values);
  for (int j = 0; j < indices.size(); j++)
    // keep diagonal
    if (indices[j] == localRowToZero)
      newvalues[j] = values[j];
    else
      newvalues[j] = Teuchos::as<Scalar>(0.25);
  A->replaceLocalValues(localRowToZero, indices, newvalues);

  // row 5 should not be Dirichlet
  drows = Utilities::DetectDirichletRows(*A, TST::magnitude(0.24));
  TEST_EQUALITY(drows[localRowToZero], false);
  TEST_EQUALITY(drows[localRowToZero - 1], false);

  // row 5 should be Dirichlet
  drows = Utilities::DetectDirichletRows(*A, TST::magnitude(0.26));
  TEST_EQUALITY(drows[localRowToZero], true);
  TEST_EQUALITY(drows[localRowToZero - 1], false);

}  // DetectDirichletRows

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, GetDiagonalInverse, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

  // blocked diagonal operator (Xpetra)
  RCP<Matrix> A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);

  RCP<const BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
  TEST_EQUALITY(bA != Teuchos::null, true);

  {
    RCP<Vector> diagInv = Utilities::GetMatrixDiagonalInverse(*bA);
    // diagInv->describe(out,Teuchos::VERB_EXTREME);
    RCP<BlockedVector> bDiagInv = Teuchos::rcp_dynamic_cast<BlockedVector>(diagInv);
    TEST_EQUALITY(bDiagInv.is_null(), false);
    TEST_EQUALITY(bDiagInv->getBlockedMap()->isSameAs(*(bA->getRangeMap())), true);
    RCP<MultiVector> diagInvMerged              = bDiagInv->Merge();
    Teuchos::ArrayRCP<const Scalar> diagInvData = diagInvMerged->getData(0);
    for (size_t i = 0; i < Teuchos::as<size_t>(diagInvData.size()); ++i) {
      if (i < 5) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0), true);
      if (i >= 5 && i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5), true);
      if (i >= 10 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0 / 3.0), true);
    }
  }

  A  = Teuchos::null;
  bA = Teuchos::null;

  // blocked diagonal operator (Thyra)
  A = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);

  bA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
  TEST_EQUALITY(bA != Teuchos::null, true);
  {
    RCP<Vector> diagInv         = Utilities::GetMatrixDiagonalInverse(*bA);
    RCP<BlockedVector> bDiagInv = Teuchos::rcp_dynamic_cast<BlockedVector>(diagInv);
    TEST_EQUALITY(bDiagInv.is_null(), false);
    TEST_EQUALITY(bDiagInv->getBlockedMap()->isSameAs(*(bA->getRangeMap())), true);
    RCP<MultiVector> diagInvMerged              = bDiagInv->Merge();
    Teuchos::ArrayRCP<const Scalar> diagInvData = diagInvMerged->getData(0);
    for (size_t i = 0; i < Teuchos::as<size_t>(diagInvData.size()); ++i) {
      if (i < 5) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0), true);
      if (i >= 5 && i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5), true);
      if (i >= 10 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0 / 3.0), true);
    }
  }
  // reordered blocked diagonal operator (Xpetra)
  A                                                   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrix(lib, 3, comm);
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm = Xpetra::blockedReorderFromString("[ [ 2 0] 1 ]");
  RCP<const BlockedCrsMatrix> bAA                     = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
  bA                                                  = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(buildReorderedBlockedCrsMatrix(brm, bAA));

  TEST_EQUALITY(bA->Rows(), 2);
  TEST_EQUALITY(bA->Cols(), 2);

  TEST_EQUALITY(bA->getRangeMapExtractor()->getThyraMode(), false);
  TEST_EQUALITY(bA->getDomainMapExtractor()->getThyraMode(), false);
  {
    RCP<Vector> diagInv         = Utilities::GetMatrixDiagonalInverse(*bA);
    RCP<BlockedVector> bDiagInv = Teuchos::rcp_dynamic_cast<BlockedVector>(diagInv);
    TEST_EQUALITY(bDiagInv.is_null(), false);
    TEST_EQUALITY(bDiagInv->getBlockedMap()->isSameAs(*(bA->getRangeMap())), true);
    RCP<MultiVector> diagInvMerged              = bDiagInv->Merge();
    Teuchos::ArrayRCP<const Scalar> diagInvData = diagInvMerged->getData(0);
    for (size_t i = 0; i < Teuchos::as<size_t>(diagInvData.size()); ++i) {
      if (i >= 10 && i < 15) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0), true);
      if (i >= 15 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5), true);
      if (i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0 / 3.0), true);
    }
  }
  // reordered blocked diagonal operator (Thyra)
  A   = TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateBlockDiagonalExampleMatrixThyra(lib, 3, comm);
  brm = Xpetra::blockedReorderFromString("[ [ 2 0] 1 ]");
  bAA = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
  bA  = Teuchos::rcp_dynamic_cast<const Xpetra::ReorderedBlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(buildReorderedBlockedCrsMatrix(brm, bAA));

  TEST_EQUALITY(bA->Rows(), 2);
  TEST_EQUALITY(bA->Cols(), 2);

  TEST_EQUALITY(bA->getRangeMapExtractor()->getThyraMode(), true);
  TEST_EQUALITY(bA->getDomainMapExtractor()->getThyraMode(), true);
  {
    RCP<Vector> diagInv         = Utilities::GetMatrixDiagonalInverse(*bA);
    RCP<BlockedVector> bDiagInv = Teuchos::rcp_dynamic_cast<BlockedVector>(diagInv);
    TEST_EQUALITY(bDiagInv.is_null(), false);
    TEST_EQUALITY(bDiagInv->getBlockedMap()->isSameAs(*(bA->getRangeMap())), true);
    RCP<MultiVector> diagInvMerged              = bDiagInv->Merge();
    Teuchos::ArrayRCP<const Scalar> diagInvData = diagInvMerged->getData(0);
    for (size_t i = 0; i < Teuchos::as<size_t>(diagInvData.size()); ++i) {
      if (i >= 10 && i < 15) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0), true);
      if (i >= 15 && i < 20) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(0.5), true);
      if (i < 10) TEST_EQUALITY(diagInvData[i] == Teuchos::as<Scalar>(1.0 / 3.0), true);
    }
  }
}  // GetDiagonalInverse

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, GetLumpedDiagonal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // Note: Lumped diagonal does not support blocked operations.
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

  RCP<CrsMatrix> crs_mat;
  RCP<const Map> row_map, col_map;

  {
    ArrayRCP<size_t> nnzPerRow(9, 0);
    row_map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9 * comm->getSize(), 0, comm);
    if (comm->getSize() == 1) {
      col_map      = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 15, 0, comm);
      nnzPerRow[0] = 3;
      nnzPerRow[1] = 4;
      nnzPerRow[2] = 4;
      nnzPerRow[3] = 4;
      nnzPerRow[4] = 5;
      nnzPerRow[5] = 5;
      nnzPerRow[6] = 4;
      nnzPerRow[7] = 5;
      nnzPerRow[8] = 5;
    } else if (comm->getSize() == 4) {
      Teuchos::Array<GlobalOrdinal> global_colids(15);
      const GlobalOrdinal rank_offset = static_cast<GlobalOrdinal>(comm->getRank() * 9);
      for (GlobalOrdinal colIdx = 0; colIdx < 9; ++colIdx) {
        global_colids[colIdx] = colIdx + rank_offset;
      }

      if (comm->getRank() == 0) {
        global_colids[9]  = 9;
        global_colids[10] = 12;
        global_colids[11] = 15;
        global_colids[12] = 18;
        global_colids[13] = 19;
        global_colids[14] = 20;

        nnzPerRow[0] = 3;
        nnzPerRow[1] = 4;
        nnzPerRow[2] = 4;
        nnzPerRow[3] = 4;
        nnzPerRow[4] = 5;
        nnzPerRow[5] = 5;
        nnzPerRow[6] = 4;
        nnzPerRow[7] = 5;
        nnzPerRow[8] = 5;
      } else if (comm->getRank() == 1) {
        global_colids[9]  = 2;
        global_colids[10] = 5;
        global_colids[11] = 8;
        global_colids[12] = 27;
        global_colids[13] = 28;
        global_colids[14] = 29;

        nnzPerRow[0] = 4;
        nnzPerRow[1] = 4;
        nnzPerRow[2] = 3;
        nnzPerRow[3] = 5;
        nnzPerRow[4] = 5;
        nnzPerRow[5] = 4;
        nnzPerRow[6] = 5;
        nnzPerRow[7] = 5;
        nnzPerRow[8] = 4;
      } else if (comm->getRank() == 2) {
        global_colids[9]  = 6;
        global_colids[10] = 7;
        global_colids[11] = 8;
        global_colids[12] = 27;
        global_colids[13] = 30;
        global_colids[14] = 33;

        nnzPerRow[0] = 4;
        nnzPerRow[1] = 5;
        nnzPerRow[2] = 5;
        nnzPerRow[3] = 4;
        nnzPerRow[4] = 5;
        nnzPerRow[5] = 5;
        nnzPerRow[6] = 3;
        nnzPerRow[7] = 4;
        nnzPerRow[8] = 4;
      } else if (comm->getRank() == 3) {
        global_colids[9]  = 15;
        global_colids[10] = 16;
        global_colids[11] = 17;
        global_colids[12] = 20;
        global_colids[13] = 23;
        global_colids[14] = 26;

        nnzPerRow[0] = 5;
        nnzPerRow[1] = 5;
        nnzPerRow[2] = 4;
        nnzPerRow[3] = 5;
        nnzPerRow[4] = 5;
        nnzPerRow[5] = 4;
        nnzPerRow[6] = 4;
        nnzPerRow[7] = 4;
        nnzPerRow[8] = 3;
      }
      col_map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 36, global_colids.view(0, 15), 0, comm);
    }
    crs_mat = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(row_map, col_map, nnzPerRow);
  }

  // Some numbers we need...
  const Scalar zero  = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one   = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two   = one + one;
  const Scalar four  = one + one + one + one;
  const Scalar six   = four + two;
  const Scalar seven = six + one;
  const Scalar eight = seven + one;
  const Scalar ten   = four + four + two;

  {
    using size_type = typename ArrayRCP<size_t>::size_type;
    ArrayRCP<size_t> row_ptr;
    ArrayRCP<LocalOrdinal> colInds;
    ArrayRCP<Scalar> values;

    if (comm->getRank() == 0) {
      size_t raw_row_ptr[] = {0, 3, 7, 11, 15, 20, 25, 29, 34, 39};
      ArrayView<size_t> row_ptr_view(raw_row_ptr, 10, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      row_ptr.deepCopy(row_ptr_view.getConst());

      LocalOrdinal raw_colInds[] = {0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 9, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 10, 3, 6, 7, 12, 4, 6, 7, 8, 13, 5, 7, 8, 10, 14};
      ArrayView<LocalOrdinal> colInds_view(raw_colInds, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      colInds.deepCopy(colInds_view.getConst());

      Scalar raw_values[] = {four, -one, -one,
                             -one, four, -one, -one,
                             -one, four, -one, -one,
                             -one, four, -one, -one,
                             -one, -one / four, four, -one / four, -one,
                             zero, zero, two, zero, zero,
                             -one / ten, four / ten, -one / ten, -one / ten,
                             -one, -one, four, -one, -one,
                             -one, -one, four, -one, -one};
      ArrayView<Scalar> values_view(raw_values, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      values.deepCopy(values_view.getConst());
    } else if (comm->getRank() == 1) {
      size_t raw_row_ptr[] = {0, 3, 7, 11, 15, 20, 25, 29, 34, 39};
      ArrayView<size_t> row_ptr_view(raw_row_ptr, 10, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      row_ptr.deepCopy(row_ptr_view.getConst());

      LocalOrdinal raw_colInds[] = {0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 9, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 10, 3, 6, 7, 12, 4, 6, 7, 8, 13, 5, 7, 8, 10, 14};
      ArrayView<LocalOrdinal> colInds_view(raw_colInds, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      colInds.deepCopy(colInds_view.getConst());

      Scalar raw_values[] = {four, -one, -one, -one, four, -one, -one, -one, four, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one};
      ArrayView<Scalar> values_view(raw_values, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      values.deepCopy(values_view.getConst());
    } else if (comm->getRank() == 2) {
      size_t raw_row_ptr[] = {0, 3, 7, 11, 15, 20, 25, 29, 34, 39};
      ArrayView<size_t> row_ptr_view(raw_row_ptr, 10, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      row_ptr.deepCopy(row_ptr_view.getConst());

      LocalOrdinal raw_colInds[] = {0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 9, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 10, 3, 6, 7, 12, 4, 6, 7, 8, 13, 5, 7, 8, 10, 14};
      ArrayView<LocalOrdinal> colInds_view(raw_colInds, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      colInds.deepCopy(colInds_view.getConst());

      Scalar raw_values[] = {four, -one, -one, -one, four, -one, -one, -one, four, -one, -one,
                             -one, four, -one, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one,
                             -one, four, -one, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one};
      ArrayView<Scalar> values_view(raw_values, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      values.deepCopy(values_view.getConst());
    } else if (comm->getRank() == 3) {
      size_t raw_row_ptr[] = {0, 3, 7, 11, 15, 20, 25, 29, 34, 39};
      ArrayView<size_t> row_ptr_view(raw_row_ptr, 10, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      row_ptr.deepCopy(row_ptr_view.getConst());

      LocalOrdinal raw_colInds[] = {0, 1, 3, 0, 1, 2, 4, 1, 2, 5, 9, 0, 3, 4, 6, 1, 3, 4, 5, 7, 2, 4, 5, 8, 10, 3, 6, 7, 12, 4, 6, 7, 8, 13, 5, 7, 8, 10, 14};
      ArrayView<LocalOrdinal> colInds_view(raw_colInds, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      colInds.deepCopy(colInds_view.getConst());

      Scalar raw_values[] = {four, -one, -one, -one, four, -one, -one, -one, four, -one, -one,
                             -one, four, -one, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one,
                             -one, four, -one, -one, -one, -one, four, -one, -one, -one, -one, four, -one, -one};
      ArrayView<Scalar> values_view(raw_values, 39, Teuchos::ERCPNodeLookup::RCP_ENABLE_NODE_LOOKUP);
      values.deepCopy(values_view.getConst());
    }
    crs_mat->setAllValues(row_ptr, colInds, values);
  }

  RCP<Matrix> mat = Teuchos::rcp(new CrsMatrixWrap(crs_mat));

  {  // Regular lumped diag
    out << std::endl
        << "== regular lumped diagonal ==" << std::endl;
    RCP<Vector> diagLumped = Utilities::GetLumpedMatrixDiagonal(*mat);
    ArrayRCP<Scalar> diag  = diagLumped->getDataNonConst(0);
    if (comm->getRank() == 0) {
      Scalar diag_ref[] = {six, seven, seven, seven, six + one / two, two, seven / ten, eight, eight};
      for (int idx = 0; idx < 9; ++idx) {
        out << std::endl
            << "comparing diag[" << idx << "] and diag_ref[" << idx << "]";
        TEST_FLOATING_EQUALITY(diag[idx], diag_ref[idx], 100 * Teuchos::ScalarTraits<Scalar>::eps());
      }
    }
  }

  {  // doReciprocal
    out << std::endl
        << "== reciprocal lumped diagonal ==" << std::endl;
    RCP<Vector> diagLumped = Utilities::GetLumpedMatrixDiagonal(*mat, true);
    ArrayRCP<Scalar> diag  = diagLumped->getDataNonConst(0);
    if (comm->getRank() == 0) {
      Scalar diag_ref[] = {one / six, one / seven, one / seven,
                           one / seven, one / (six + one / two), one / four,
                           ten / seven, one / eight, one / eight};
      for (int idx = 0; idx < 9; ++idx) {
        out << std::endl
            << "comparing diag[" << idx << "] and diag_ref[" << idx << "]";
        TEST_FLOATING_EQUALITY(diag[idx], diag_ref[idx], 100 * Teuchos::ScalarTraits<Scalar>::eps());
      }
    }
  }

  {  // val < 0.9 treated as 0.0
    out << std::endl
        << "== reciprocal lumped diagonal with dropping ==" << std::endl;
    RCP<Vector> diagLumped = Utilities::GetLumpedMatrixDiagonal(*mat, true, 0.9);
    ArrayRCP<Scalar> diag  = diagLumped->getDataNonConst(0);
    if (comm->getRank() == 0) {
      Scalar diag_ref[] = {one / six, one / seven, one / seven,
                           one / seven, one / (six + one / two), one / four,
                           zero, one / eight, one / eight};
      for (int idx = 0; idx < 9; ++idx) {
        out << std::endl
            << "comparing diag[" << idx << "] and diag_ref[" << idx << "]";
        TEST_FLOATING_EQUALITY(diag[idx], diag_ref[idx], 100 * Teuchos::ScalarTraits<Scalar>::eps());
      }
    }
  }

  {  // nnzPerRow(i) <= 1 --> diag(i) = zero
    out << std::endl
        << "== reciprocal lumped diagonal with single-entry row ==" << std::endl;
    RCP<Vector> diagLumped = Utilities::GetLumpedMatrixDiagonal(*mat, true, 0, 42, true);
    ArrayRCP<Scalar> diag  = diagLumped->getDataNonConst(0);
    if (comm->getRank() == 0) {
      Scalar diag_ref[] = {one / six, one / seven, one / seven,
                           one / seven, one / (six + one / two), zero,
                           ten / seven, one / eight, one / eight};
      for (int idx = 0; idx < 9; ++idx) {
        out << std::endl
            << "comparing diag[" << idx << "] and diag_ref[" << idx << "]";
        TEST_FLOATING_EQUALITY(diag[idx], diag_ref[idx], 100 * Teuchos::ScalarTraits<Scalar>::eps());
      }
    }
  }

  {  // nnzPerRow(i) <= 1 --> diag(i) = zero
    out << std::endl
        << "== reciprocal lumped diagonal with zero-entry row ==" << std::endl;
    RCP<Vector> diagLumped = Utilities::GetLumpedMatrixDiagonal(*mat, true, 0.9, 42, true);
    ArrayRCP<Scalar> diag  = diagLumped->getDataNonConst(0);
    if (comm->getRank() == 0) {
      Scalar diag_ref[] = {one / six, one / seven, one / seven,
                           one / seven, one / (six + one / two), zero,
                           ten + ten + ten + ten + two, one / eight, one / eight};
      for (int idx = 0; idx < 9; ++idx) {
        out << std::endl
            << "comparing diag[" << idx << "] and diag_ref[" << idx << "]";
        TEST_FLOATING_EQUALITY(diag[idx], diag_ref[idx], 100 * Teuchos::ScalarTraits<Scalar>::eps());
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, GetInverse, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar> TST;

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();

  Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

  RCP<Map> m = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 100, 0, comm);

  RCP<Vector> v  = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(m, true);
  RCP<Vector> tv = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(m, true);
  {
    Teuchos::ArrayRCP<Scalar> vData  = v->getDataNonConst(0);
    Teuchos::ArrayRCP<Scalar> tvData = tv->getDataNonConst(0);
    for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(v->getLocalLength()); ++i) {
      vData[i]  = Teuchos::as<Scalar>(i + 1);
      tvData[i] = Teuchos::ScalarTraits<Scalar>::one() / Teuchos::as<Scalar>(i + 1);
    }
  }
  RCP<Vector> inv = Utilities::GetInverse(v);

  tv->update(1.0, *inv, -1.0);
  TEST_EQUALITY(tv->norm1(), Teuchos::ScalarTraits<Scalar>::zero());
  TEST_EQUALITY(tv->norm2(), Teuchos::ScalarTraits<Scalar>::zero());
  TEST_EQUALITY(tv->normInf(), Teuchos::ScalarTraits<Scalar>::zero());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, GetThresholdedMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST            = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename TST::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();

  // Don't test for complex - matrix reader won't work
  if (TST::isComplex) {
    success = true;
    return;
  }
  RCP<Matrix> Ain = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("TestMatrices/filter.mm", lib, comm);

  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Aout =
      MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetThresholdedMatrix(Ain, 1e-5, true, -1);

  TEST_EQUALITY(Aout->getCrsGraph()->getGlobalNumEntries(), Teuchos::as<size_t>(13));
  TEST_FLOATING_EQUALITY(Aout->getFrobeniusNorm(), Teuchos::as<magnitude_type>(7.549834435270750), 1e2 * Teuchos::ScalarTraits<Scalar>::eps());
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, GetThresholdedGraph, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  using TST            = Teuchos::ScalarTraits<SC>;
  using magnitude_type = typename TST::magnitudeType;
  using TMT            = Teuchos::ScalarTraits<magnitude_type>;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = MueLuTests::TestHelpers::Parameters::getLib();

  // Don't test for complex - matrix reader won't work
  if (TST::isComplex) {
    success = true;
    return;
  }
  RCP<Matrix> A = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("TestMatrices/filter.mm", lib, comm);

  RCP<Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> graph =
      MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetThresholdedGraph(A, 1e-5, -1);

  TEST_EQUALITY(graph->getGlobalNumEntries(), Teuchos::as<size_t>(13));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Utilities, TransposeNonsymmetricConstMatrix, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  RCP<const Teuchos::Comm<int>> comm = TestHelpers::Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers::Parameters::getLib();

  const GO numGlobalElements = 29;
  RCP<Map> dofMap            = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, numGlobalElements, 0, comm);

  TEST_ASSERT(!dofMap.is_null());
  TEST_EQUALITY_CONST(dofMap->getGlobalNumElements(), numGlobalElements);

  RCP<const Matrix> matrix = TestHelpers::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 1.0, 2.0, 3.0, lib);
  TEST_ASSERT(!matrix.is_null());
  TEST_EQUALITY_CONST(matrix->getGlobalNumRows(), numGlobalElements);

  RCP<const Matrix> transposedMatrix = Utilities::Transpose(const_cast<Matrix&>(*matrix));
  TEST_ASSERT(!transposedMatrix.is_null());

  TEST_ASSERT(transposedMatrix->getRangeMap()->isSameAs(*matrix->getDomainMap()));
  TEST_ASSERT(transposedMatrix->getDomainMap()->isSameAs(*matrix->getRangeMap()));

  // Verify, that A^T actually differs from A
  {
    RCP<Matrix> diffMatrix = rcp(new CrsMatrixWrap(matrix->getCrsGraph()));
    MatrixMatrix::TwoMatrixAdd(*matrix, false, 1.0, *transposedMatrix, false, -1.0, diffMatrix, out);
    diffMatrix->fillComplete();
    TEST_ASSERT(!diffMatrix.is_null());

    bool allEntriesAreZero = true;
    for (LO lRowId = 0; lRowId < Teuchos::as<LO>(diffMatrix->getLocalNumRows()); ++lRowId) {
      ArrayView<const LO> cols;
      ArrayView<const Scalar> vals;
      diffMatrix->getLocalRowView(lRowId, cols, vals);

      TEST_INEQUALITY_CONST(cols.size(), 0);
      TEST_INEQUALITY_CONST(vals.size(), 0);

      for (const auto& entry : vals) {
        if (entry != Teuchos::ScalarTraits<Scalar>::zero()) allEntriesAreZero = false;
      }
    }
    TEST_ASSERT(!allEntriesAreZero);
  }

  RCP<const Matrix> doubleTransposedMatrix = Utilities::Transpose(const_cast<Matrix&>(*transposedMatrix));
  TEST_ASSERT(!doubleTransposedMatrix.is_null());

  TEST_ASSERT(doubleTransposedMatrix->getRangeMap()->isSameAs(*matrix->getRangeMap()));
  TEST_ASSERT(doubleTransposedMatrix->getDomainMap()->isSameAs(*matrix->getDomainMap()));

  // Transpose twice: A - (A^T)^T needs to be the zero matrix
  {
    RCP<Matrix> diffMatrix = rcp(new CrsMatrixWrap(matrix->getCrsGraph()));
    MatrixMatrix::TwoMatrixAdd(*matrix, false, 1.0, *doubleTransposedMatrix, false, -1.0, diffMatrix, out);
    diffMatrix->fillComplete();
    TEST_ASSERT(!diffMatrix.is_null());

    bool allEntriesAreZero = true;
    for (LO lRowId = 0; lRowId < Teuchos::as<LO>(diffMatrix->getLocalNumRows()); ++lRowId) {
      ArrayView<const LO> cols;
      ArrayView<const Scalar> vals;
      diffMatrix->getLocalRowView(lRowId, cols, vals);

      TEST_INEQUALITY_CONST(cols.size(), 0);
      TEST_INEQUALITY_CONST(vals.size(), 0);

      for (const auto& entry : vals) {
        if (entry != Teuchos::ScalarTraits<Scalar>::zero()) allEntriesAreZero = false;
      }
    }
    TEST_ASSERT(allEntriesAreZero);
  }
}

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, MatMatMult_EpetraVsTpetra, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, DetectDirichletRows, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, GetDiagonalInverse, Scalar, LO, GO, Node)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, GetLumpedDiagonal, Scalar, LO, GO, Node)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, GetInverse, Scalar, LO, GO, Node)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, GetThresholdedMatrix, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, GetThresholdedGraph, Scalar, LO, GO, Node)       \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Utilities, TransposeNonsymmetricConstMatrix, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
