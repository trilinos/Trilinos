// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "SingletonFiltering_TestUtils.hpp"

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHelpers.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>

namespace {

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Tpetra::TestingUtilities::getDefaultComm;

// Test Class for SF1 (SingletonFiltering Problem 1)
// --------------------------------------------------------------------------
// Derived test class so it is easy to access protected data and functions.
template <typename Scalar, typename LO, typename GO, typename Node>
class Test_CrsSingletonFilter_LinearProblem
  : public Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> {
 public:
  Test_CrsSingletonFilter_LinearProblem(bool run_on_host = false, bool verbose = false)
    : Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>(run_on_host, verbose) {}

  void setFullMatrix(typename Teuchos::RCP<typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::row_matrix_type> matrix) { this->FullMatrix_ = matrix; }

  void test_InitFullMatrixAccess() {
    this->InitFullMatrixAccess();

    TEUCHOS_ASSERT(this->localMaxNumRowEntries_ == 6);
    TEUCHOS_ASSERT(this->FullMatrixIsCrsMatrix_ == true);
    TEUCHOS_ASSERT(this->FullCrsMatrix_ != Teuchos::null);
  }

  void test_GetRow3(Teuchos::FancyOStream &out, bool &success) {  // 3 arguments
    using LocalOrdinal_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::local_ordinal_type;
    size_t NumIndices    = 1;
    Teuchos::Array<LocalOrdinal_t> localIndices;
    // Row 0
    this->GetRow(0, NumIndices, localIndices);
    Teuchos::Array<LocalOrdinal_t> ansIndices = Teuchos::Array<LocalOrdinal_t>({1, 2});
    TEUCHOS_ASSERT(NumIndices == static_cast<size_t>(ansIndices.size()));
    TEST_COMPARE_ARRAYS(localIndices, ansIndices);
    // Row 8
    this->GetRow(8, NumIndices, localIndices);
    ansIndices = Teuchos::Array<LocalOrdinal_t>({6, 8, 9, 10, 12, 13});
    TEUCHOS_ASSERT(NumIndices == static_cast<size_t>(ansIndices.size()));
    TEST_COMPARE_ARRAYS(localIndices, ansIndices);
    // Row 10 - singleton
    this->GetRow(10, NumIndices, localIndices);
    ansIndices = Teuchos::Array<LocalOrdinal_t>({8});
    TEUCHOS_ASSERT(NumIndices == static_cast<size_t>(ansIndices.size()));
    TEST_COMPARE_ARRAYS(localIndices, ansIndices);
  }

  void test_GetRow4(Teuchos::FancyOStream &out, bool &success) {  // 4 arguments
    using LocalOrdinal_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::local_ordinal_type;
    size_t NumEntries    = 0;
    Teuchos::ArrayView<const Scalar> Values;
    Teuchos::ArrayView<const LocalOrdinal_t> localIndices;
    // Row 8
    this->GetRow(8, NumEntries, Values, localIndices);
    Teuchos::Array<LocalOrdinal_t> ansIndices = Teuchos::Array<LocalOrdinal_t>({6, 8, 9, 10, 12, 13});
    TEST_COMPARE_ARRAYS(localIndices, ansIndices);
    double relTol                    = 1.0e-05;
    Teuchos::Array<Scalar> ansValues = Teuchos::Array<Scalar>({Scalar(-1.0), Scalar(0.00140635), Scalar(-1.0), Scalar(1.0), Scalar(-0.0012461), Scalar(-0.000160258)});
    if constexpr (std::is_same<Scalar, long long>::value) {
      // It appears that Reader_t::readSparseFile() and Scalar() round/truncate differently for long long,
      // so need to explicitly set ansValues to match Reader_t::readSparseFile().
      ansValues = Teuchos::Array<Scalar>({Scalar(-1), Scalar(1), Scalar(-1), Scalar(1), Scalar(-1), Scalar(-1)});
    }
    TEST_COMPARE_FLOATING_ARRAYS(Values, ansValues, relTol);
  }

  void test_GetRowGCIDs(Teuchos::FancyOStream &out, bool &success) {
    using GlobalOrdinal_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::global_ordinal_type;
    size_t NumEntries     = 0;
    Teuchos::ArrayView<const Scalar> Values;
    Teuchos::Array<GlobalOrdinal_t> Indices;
    // Row 8
    this->GetRowGCIDs(8, NumEntries, Values, Indices);
    Teuchos::Array<GlobalOrdinal_t> ansIndices = Teuchos::Array<GlobalOrdinal_t>({6, 8, 9, 10, 12, 13});
    TEST_COMPARE_ARRAYS(Indices, ansIndices);
    double relTol                    = 1.0e-05;
    Teuchos::Array<Scalar> ansValues = Teuchos::Array<Scalar>({Scalar(-1.0), Scalar(0.00140635), Scalar(-1.0), Scalar(1.0), Scalar(-0.0012461), Scalar(-0.000160258)});
    if constexpr (std::is_same<Scalar, long long>::value) {
      // It appears that Reader_t::readSparseFile() and Scalar() round/truncate differently for long long,
      // so need to explicitly set ansValues to match Reader_t::readSparseFile().
      ansValues = Teuchos::Array<Scalar>({Scalar(-1), Scalar(1), Scalar(-1), Scalar(1), Scalar(-1), Scalar(-1)});
    }
    TEST_COMPARE_FLOATING_ARRAYS(Values, ansValues, relTol);
  }

  void test_GenerateReducedMap(Teuchos::FancyOStream &out, bool &success) {
    using Map_t       = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::map_type;
    using VectorInt_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::vector_type_int;

    Teuchos::RCP<VectorInt_t> mapColors = rcp(new VectorInt_t(this->FullMatrixRowMap()));
    Teuchos::Array<int> colors          = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < static_cast<size_t>(colors.size()); ++i) {
      mapColors->replaceGlobalValue(static_cast<int>(i), colors[i]);
    }

    Teuchos::RCP<const Map_t> reducedMap = this->GenerateReducedMap(this->FullMatrixRowMap(), mapColors, 0);
    auto Comm                            = Tpetra::getDefaultComm();
    Teuchos::Array<GO> globalIndices     = {0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 15, 16, 17, 18};
    GO numGlobalElements                 = 17;
    GO indexBase                         = 0;
    Teuchos::RCP<Map_t> ansReducedMap    = Teuchos::rcp(new Map_t(numGlobalElements, globalIndices(), indexBase, Comm));

    TEST_ASSERT(reducedMap->isSameAs(*ansReducedMap));
  }

  void test_Analyze(Teuchos::FancyOStream &out, bool &success) {
    using LocalOrdinal_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::local_ordinal_type;
    using VectorInt_t    = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::vector_type_int;

    TEST_ASSERT(this->origObj_ != Teuchos::null);
    TEST_ASSERT(this->FullProblem() == Teuchos::null);
    TEST_ASSERT(this->AnalysisDone_ == true);
    TEST_ASSERT(this->FullMatrix()->getGlobalNumRows() == 19);
    TEST_ASSERT(this->FullMatrix()->getGlobalNumEntries() == 54);

    // Check color maps
    Teuchos::RCP<VectorInt_t> mapColors = rcp(new VectorInt_t(this->FullMatrixRowMap()));
    Teuchos::Array<int> colors          = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < static_cast<size_t>(colors.size()); ++i) {
      mapColors->replaceGlobalValue(static_cast<int>(i), colors[i]);
    }
    auto rowMapColorsData = this->RowMapColors_->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto mapColorsData    = mapColors->getLocalViewHost(Tpetra::Access::ReadOnly);
    TEST_ASSERT(rowMapColorsData.extent(0) == mapColorsData.extent(0));
    TEST_ASSERT(rowMapColorsData.extent(1) == mapColorsData.extent(1));
    TEST_ASSERT(rowMapColorsData.extent(1) == 1);
    for (size_t i = 0; i < rowMapColorsData.extent(0); ++i) {
      TEST_ASSERT(rowMapColorsData(i, 0) == mapColorsData(i, 0));
    }
    auto colMapColorsData = this->ColMapColors_->getLocalViewHost(Tpetra::Access::ReadOnly);
    TEST_ASSERT(colMapColorsData.extent(0) == mapColorsData.extent(0));
    TEST_ASSERT(colMapColorsData.extent(1) == mapColorsData.extent(1));
    TEST_ASSERT(colMapColorsData.extent(1) == 1);
    for (size_t i = 0; i < colMapColorsData.extent(0); ++i) {
      TEST_ASSERT(colMapColorsData(i, 0) == mapColorsData(i, 0));
    }

    TEST_ASSERT(this->localNumSingletonRows_ == 1);
    TEST_ASSERT(this->localNumSingletonCols_ == 1);

    {
      auto h_ColSingletonRowLIDs = Kokkos::create_mirror_view(this->ColSingletonRowLIDs_);
      Kokkos::deep_copy(h_ColSingletonRowLIDs, this->ColSingletonRowLIDs_);
      Teuchos::ArrayRCP<LocalOrdinal_t> ColSingletonRowLIDs    = Teuchos::arcp(h_ColSingletonRowLIDs.data(), 0, h_ColSingletonRowLIDs.extent(0), false);
      Teuchos::ArrayRCP<LocalOrdinal_t> ansColSingletonRowLIDs = Teuchos::arcp(new LocalOrdinal_t[1]{8}, 0, 1, true);
      TEST_COMPARE_ARRAYS(ColSingletonRowLIDs, ansColSingletonRowLIDs);
    }
    {
      auto h_ColSingletonColLIDs = Kokkos::create_mirror_view(this->ColSingletonColLIDs_);
      Kokkos::deep_copy(h_ColSingletonColLIDs, this->ColSingletonColLIDs_);
      Teuchos::ArrayRCP<LocalOrdinal_t> ColSingletonColLIDs    = Teuchos::arcp(h_ColSingletonColLIDs.data(), 0, h_ColSingletonColLIDs.extent(0), false);
      Teuchos::ArrayRCP<LocalOrdinal_t> ansColSingletonColLIDs = Teuchos::arcp(new LocalOrdinal_t[1]{10}, 0, 1, true);
      TEST_COMPARE_ARRAYS(ColSingletonColLIDs, ansColSingletonColLIDs);
    }
  }

  void test_Operator(Teuchos::FancyOStream &out, bool &success) {
    using CrsMatrix_t   = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::crs_matrix_type;
    using MultiVector_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::multivector_type;
    using Reader_t      = typename Tpetra::MatrixMarket::Reader<CrsMatrix_t>;

    // Tests for ConstructReducedProblem()
    TEST_ASSERT(this->HaveReducedProblem_ == true);
    TEST_ASSERT(this->FullProblem() != Teuchos::null);
    TEST_ASSERT(this->FullMatrix() != Teuchos::null);
    TEST_ASSERT(this->FullProblem()->getRHS() != Teuchos::null);
    TEST_ASSERT(this->FullProblem()->getLHS() != Teuchos::null);
    TEST_ASSERT(this->SingletonsDetected() == true);

    auto Comm = Tpetra::getDefaultComm();

    auto reducedRowMap    = this->ReducedMatrix()->getRowMap();
    auto reducedColMap    = this->ReducedMatrix()->getColMap();
    auto reducedDomainMap = this->ReducedMatrix()->getDomainMap();
    auto reducedRangeMap  = this->ReducedMatrix()->getRangeMap();

    RCP<CrsMatrix_t> A_Reduced;
    RCP<MultiVector_t> LHS_Reduced, RHS_Reduced;
    A_Reduced   = Reader_t::readSparseFile("SF1_Matrix_Reduced.mm", reducedRowMap, reducedColMap, reducedDomainMap, reducedRangeMap);
    LHS_Reduced = Reader_t::readDenseFile("SF1_LHS_Reduced.mm", Comm, reducedRowMap);
    RHS_Reduced = Reader_t::readDenseFile("SF1_RHS_Reduced.mm", Comm, reducedRowMap);

    TEUCHOS_ASSERT(compareCrsMatrices(
        Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(A_Reduced, true),
        Teuchos::rcp_dynamic_cast<const CrsMatrix_t>(this->ReducedMatrix(), true), Comm, out, 1.0e-05));

    TEUCHOS_ASSERT(compareMultiVectors(
        Teuchos::rcp_dynamic_cast<const MultiVector_t>(LHS_Reduced, true),
        Teuchos::rcp_dynamic_cast<const MultiVector_t>(this->ReducedLHS_, true), Comm, out));

    TEUCHOS_ASSERT(compareMultiVectors(
        Teuchos::rcp_dynamic_cast<const MultiVector_t>(RHS_Reduced, true),
        Teuchos::rcp_dynamic_cast<const MultiVector_t>(this->ReducedRHS_, true), Comm, out));

    // Tests for ReducedProblem()
    TEST_ASSERT(this->newObj_ != Teuchos::null);

    TEST_ASSERT(this->ReducedMatrix()->getGlobalNumRows() == A_Reduced->getGlobalNumRows());
    TEST_ASSERT(this->ReducedMatrix()->getGlobalNumEntries() == A_Reduced->getGlobalNumEntries());
  }
};

// Unit Tests
// --------------------------------------------------------------------------

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SF1, Functions, LO, GO, Scalar, Node) {
  if (Teuchos::ScalarTraits<Scalar>::isComplex) return;  // MatrixMarket reader does not work properly with complex.

  using CrsMatrix_t   = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::crs_matrix_type;
  using MultiVector_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::multivector_type;
  using Reader_t      = typename Tpetra::MatrixMarket::Reader<CrsMatrix_t>;

  auto Comm = Tpetra::getDefaultComm();

  RCP<CrsMatrix_t> A_Original;
  RCP<MultiVector_t> LHS_Original, RHS_Original;

  Test_CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> test_SF;
  A_Original = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
  test_SF.setFullMatrix(A_Original);

  test_SF.test_InitFullMatrixAccess();
  test_SF.test_GetRow3(out, success);
  test_SF.test_GetRow4(out, success);
  test_SF.test_GetRowGCIDs(out, success);
  test_SF.test_GenerateReducedMap(out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SF1, Analyze, LO, GO, Scalar, Node) {
  if (Teuchos::ScalarTraits<Scalar>::isComplex) return;  // MatrixMarket reader does not work properly with complex.

  using CrsMatrix_t     = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::crs_matrix_type;
  using MultiVector_t   = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::multivector_type;
  using Map_t           = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::map_type;
  using LinearProblem_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::linear_problem_type;
  using Reader_t        = typename Tpetra::MatrixMarket::Reader<CrsMatrix_t>;

  auto Comm = Tpetra::getDefaultComm();

  RCP<CrsMatrix_t> A_Original;
  RCP<MultiVector_t> LHS_Original, RHS_Original;

  A_Original                      = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
  RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
  LHS_Original                    = Reader_t::readDenseFile("SF1_LHS_Original.mm", Comm, A_Original_Map);
  RHS_Original                    = Reader_t::readDenseFile("SF1_RHS_Original.mm", Comm, A_Original_Map);

  RCP<MultiVector_t> x                     = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors()));
  RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t(A_Original, x, RHS_Original));

  Test_CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> test_SF;
  test_SF.analyze(preSingletonProblem);
  test_SF.test_Analyze(out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SF1, Operator, LO, GO, Scalar, Node) {
  // operator() just calls analyze(LinearProblem) and construct(), so this test will
  // basically cover construct() since analyze() is covered above.

  if (Teuchos::ScalarTraits<Scalar>::isComplex) return;  // MatrixMarket reader does not work properly with complex.

  using CrsMatrix_t     = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::crs_matrix_type;
  using MultiVector_t   = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::multivector_type;
  using Map_t           = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::map_type;
  using LinearProblem_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>::linear_problem_type;
  using Reader_t        = typename Tpetra::MatrixMarket::Reader<CrsMatrix_t>;

  auto Comm = Tpetra::getDefaultComm();

  RCP<CrsMatrix_t> A_Original;
  RCP<MultiVector_t> LHS_Original, RHS_Original;

  A_Original                      = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
  RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
  LHS_Original                    = Reader_t::readDenseFile("SF1_LHS_Original.mm", Comm, A_Original_Map);
  RHS_Original                    = Reader_t::readDenseFile("SF1_RHS_Original.mm", Comm, A_Original_Map);

  RCP<MultiVector_t> x                     = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors()));
  RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t(A_Original, x, RHS_Original));

  Test_CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node> test_SF;
  RCP<LinearProblem_t> postSingletonProblem = test_SF(preSingletonProblem);
  test_SF.test_Operator(out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SF2, fwd, LO, GO, Scalar, Node) {
  if (Teuchos::ScalarTraits<Scalar>::isComplex) return;  // MatrixMarket reader does not work properly with complex.

  auto Comm = Tpetra::getDefaultComm();

  test_Singleton_fwd<Scalar, LO, GO, Node>(
      "SF2_Matrix_Original.mm", "SF2_LHS_Original.mm", "SF2_RHS_Original.mm",
      "SF2_Matrix_Reduced.mm", "SF2_LHS_Reduced.mm", "SF2_RHS_Reduced.mm",
      Comm, out, success);
}

#define UNIT_TEST_GROUP(SCALAR, LO, GO, NODE)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SF1, Functions, LO, GO, SCALAR, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SF1, Analyze, LO, GO, SCALAR, NODE)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SF1, Operator, LO, GO, SCALAR, NODE)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SF2, fwd, LO, GO, SCALAR, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP)

}  // namespace
