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
#include <Tpetra_CrsSingletonFilter_LinearProblem.hpp>
#include <Tpetra_SolverMap_LinearProblem.hpp>
#include <Tpetra_Reindex_LinearProblem.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <TpetraExt_MatrixMatrix.hpp>

namespace {

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Tpetra::TestingUtilities::getDefaultComm;

// Unit Tests
// --------------------------------------------------------------------------

// Singleton Filtering Test
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SingletonFiltering_P1, fwd, LO, GO, Scalar, Node) {
  if (Teuchos::ScalarTraits<Scalar>::isComplex) return;  // MatrixMarket reader does not work properly with complex.

  auto Comm = Tpetra::getDefaultComm();

  test_Singleton_fwd<Scalar, LO, GO, Node>(
      "SF1_Matrix_Original.mm", "SF1_LHS_Original.mm", "SF1_RHS_Original.mm",
      "SF1_Matrix_Reduced.mm", "SF1_LHS_Reduced.mm", "SF1_RHS_Reduced.mm",
      Comm, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SingletonFiltering_P1, MultipleTransforms, LO, GO, Scalar, Node) {
  // MultipleTransforms sequences several transforms together (SingletonFiltering, SolverMap and ReIndex)

  if (Teuchos::ScalarTraits<Scalar>::isComplex) return;  // MatrixMarket reader does not work properly with complex.

  using Problem_t               = typename Tpetra::LinearProblem<Scalar, LO, GO, Node>;
  using CrsSingletonFiltering_t = typename Tpetra::CrsSingletonFilter_LinearProblem<Scalar, LO, GO, Node>;
  using CrsMatrix_t             = typename CrsSingletonFiltering_t::crs_matrix_type;
  using MultiVector_t           = typename CrsSingletonFiltering_t::multivector_type;
  using Map_t                   = typename CrsSingletonFiltering_t::map_type;
  using LinearProblem_t         = typename CrsSingletonFiltering_t::linear_problem_type;
  using SolverMap_t             = typename Tpetra::SolverMap_LinearProblem<Scalar, LO, GO, Node>;
  using Reindex_t               = typename Tpetra::Reindex_LinearProblem<Scalar, LO, GO, Node>;
  using Reader_t                = typename Tpetra::MatrixMarket::Reader<CrsMatrix_t>;

  auto Comm = Tpetra::getDefaultComm();

  RCP<CrsMatrix_t> A_Original;
  RCP<MultiVector_t> LHS_Original, RHS_Original;

  A_Original                      = Reader_t::readSparseFile("SF1_Matrix_Original.mm", Comm);
  RCP<const Map_t> A_Original_Map = A_Original->getRangeMap();
  LHS_Original                    = Reader_t::readDenseFile("SF1_LHS_Original.mm", Comm, A_Original_Map);
  RHS_Original                    = Reader_t::readDenseFile("SF1_RHS_Original.mm", Comm, A_Original_Map);

  bool verbose                             = true;
  bool run_on_host                         = false;
  RCP<MultiVector_t> x                     = rcp(new MultiVector_t(A_Original_Map, LHS_Original->getNumVectors()));
  RCP<LinearProblem_t> preSingletonProblem = rcp(new LinearProblem_t(A_Original, x, RHS_Original));

  // Singleton Transform Forward
  CrsSingletonFiltering_t SingletonTransform(run_on_host, verbose);
  RCP<Problem_t> postSingletonProblem = SingletonTransform(preSingletonProblem);
  SingletonTransform.fwd();

  // Solver Map Transform Forward
  SolverMap_t SolverMapTransform;
  RCP<Problem_t> postSolverMapProblem = SolverMapTransform(postSingletonProblem);
  SolverMapTransform.fwd();

  // Reindex Transform Forward
  Teuchos::RCP<Map_t const> newRowMap(Teuchos::null);
  Reindex_t ReindexTransform(newRowMap);
  RCP<Problem_t> postReindexProblem = ReindexTransform(postSolverMapProblem);
  ReindexTransform.fwd();

  // In lieu of a solve, just set the reduced solution to test transforms.
  auto reducedRowMap                 = SingletonTransform.ReducedMatrix()->getRowMap();
  RCP<MultiVector_t> reducedSolution = Reader_t::readDenseFile("SF1_Solution_Reduced.mm", Comm, reducedRowMap);

  postReindexProblem->getLHS()->update(1.0, *reducedSolution, 0.0);

  // Reindex Transform Reverse
  ReindexTransform.rvs();

  // Solver Map Transform Reverse
  SolverMapTransform.rvs();

  // Singleton Transform Reverse
  SingletonTransform.rvs();

  // Test so long as not long long.
  if constexpr (!std::is_same<Scalar, long long>::value) {
    RCP<MultiVector_t> solution = Reader_t::readDenseFile("SF1_Solution.mm", Comm, A_Original_Map);

    TEUCHOS_ASSERT(compareMultiVectors(
        Teuchos::rcp_dynamic_cast<const MultiVector_t>(solution, true),
        Teuchos::rcp_dynamic_cast<const MultiVector_t>(preSingletonProblem->getLHS(), true),
        Comm, out));

    out << "MultiVector solution:" << std::endl;
    solution->describe(out, Teuchos::VERB_EXTREME);
    out << "MultiVector preSingletonProblem->getLHS():" << std::endl;
    preSingletonProblem->getLHS()->describe(out, Teuchos::VERB_EXTREME);
  }
}

#define UNIT_TEST_GROUP(SCALAR, LO, GO, NODE)                                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFiltering_P1, fwd, LO, GO, SCALAR, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SingletonFiltering_P1, MultipleTransforms, LO, GO, SCALAR, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP)

}  // namespace
