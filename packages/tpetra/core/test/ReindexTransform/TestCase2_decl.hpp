// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_2_DECL_HPP
#define TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_2_DECL_HPP

#include "TestCaseBase_decl.hpp"

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
class TestCase2 : public TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> {
 public:
  using Map_t         = typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::Map_t;
  using MultiVector_t = typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::MultiVector_t;
  using Matrix_t      = typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::Matrix_t;
  using Problem_t     = typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::Problem_t;
  using ColValPair_t  = std::pair<LocalId_t /*localColIndex*/, Scalar_t>;

  TestCase2() = delete;

  TestCase2(Teuchos::RCP<Teuchos::Comm<int> const> comm);

  TestCase2(TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t> const& source) = delete;

  TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t>& operator=(TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t> const& rhs) = delete;

  ~TestCase2();

  bool checkTransformedProblem(Problem_t const* transformedProblem) const;

  bool checkAfterFwdRvs(Problem_t const* originalProblem) const;

 protected:
  void prepareDataToCreateMatrix(std::vector<std::vector<GlobalId_t> >& globalRowIndicesPerRank, std::vector<std::vector<GlobalId_t> >& globalColIndicesPerRank, std::vector<std::vector<GlobalId_t> >& globalDomainIndicesPerRank, std::vector<std::vector<GlobalId_t> >& globalRangeIndicesPerRank, std::vector<std::vector<std::vector<ColValPair_t> > >& colValuePairsPerRankAndLocalRow);
};

#include "TestCase2_def.hpp"

#endif  // TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_2_DECL_HPP
