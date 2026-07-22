// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_BASE_DECL_HPP
#define TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_BASE_DECL_HPP

#include <Tpetra_LinearProblem.hpp>

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
class TestCaseBase {
 public:
  using global_size_t = Tpetra::global_size_t;
  using Map_t         = Tpetra::Map<LocalId_t, GlobalId_t, Node_t>;
  using MultiVector_t = Tpetra::MultiVector<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Matrix_t      = Tpetra::CrsMatrix<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Problem_t     = Tpetra::LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;

  TestCaseBase() = delete;

  TestCaseBase(Teuchos::RCP<Teuchos::Comm<int> const> comm, global_size_t const globalNumRows, global_size_t const globalNumCols, global_size_t const rowIndexBase, global_size_t const colIndexBase, global_size_t const maxNnzPerRow, global_size_t const globalNumNnz, Scalar_t const frobNorm, std::vector<Scalar_t> const& lhsNorms2, std::vector<Scalar_t> const& rhsNorms2);

  TestCaseBase(TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> const& source) = delete;

  TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>& operator=(TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> const& rhs) = delete;

  virtual ~TestCaseBase();

  Problem_t* linearProblem() const;

  virtual bool checkTransformedProblem(Problem_t const* transformedProblem) const = 0;

 protected:
  void baseInstantiateLinearProblem(std::string const& matInputFileName, std::string const& lhsInputFileName, std::string const& rhsInputFileName);

  bool baseCheckBeforeOrAfterTransform(Problem_t const* problem, std::string const& caseString) const;

  bool baseCheckTransformedProblem(Problem_t const* transformedProblem) const;

  bool baseCheckAfterFwdRvs(Problem_t const* originalProblem) const;

  bool baseCheckMapsAreDifferent(Matrix_t const& mat1, Matrix_t const& mat2) const;

  bool baseCheckMatricesAreEqual(Matrix_t const& mat1, Matrix_t const& mat2) const;

  Teuchos::RCP<const Teuchos::Comm<int> > m_comm;
  int const m_numRanks;
  int const m_myRank;

  global_size_t const m_globalNumRows;
  global_size_t const m_globalNumCols;

  global_size_t const m_rowIndexBase;
  global_size_t const m_colIndexBase;

  global_size_t const m_maxNnzPerRow;
  global_size_t const m_globalNumNnz;

  Scalar_t const m_frobNorm;
  std::vector<Scalar_t> const m_lhsNorms2;
  std::vector<Scalar_t> const m_rhsNorms2;

  std::unique_ptr<Map_t> m_rowMap;
  std::unique_ptr<Map_t> m_colMap;
  std::unique_ptr<Map_t> m_domainMap;
  std::unique_ptr<Map_t> m_rangeMap;
  std::unique_ptr<Matrix_t> m_matrix;
  std::unique_ptr<MultiVector_t> m_lhs;
  std::unique_ptr<MultiVector_t> m_rhs;
  Teuchos::RCP<Matrix_t> m_matrixRCP;
  Teuchos::RCP<MultiVector_t> m_lhsRCP;
  Teuchos::RCP<MultiVector_t> m_rhsRCP;
  std::unique_ptr<Problem_t> m_linearProblem;
  bool m_mapsShallChange;
};

#include "TestCaseBase_def.hpp"

#endif  // TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_BASE_DECL_HPP
