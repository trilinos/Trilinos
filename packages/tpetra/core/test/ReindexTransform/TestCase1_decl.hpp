// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_1_DECL_HPP
#define TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_1_DECL_HPP

#include "TestCaseBase_decl.hpp"

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
class TestCase1 : public TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> {
 public:
  using Problem_t = typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::Problem_t;

  TestCase1() = delete;

  TestCase1(Teuchos::RCP<Teuchos::Comm<int> const> comm);

  TestCase1(TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t> const& source) = delete;

  TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t>& operator=(TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t> const& rhs) = delete;

  ~TestCase1();

  bool checkTransformedProblem(Problem_t const* transformedProblem) const;

  bool checkAfterFwdRvs(Problem_t const* originalProblem) const;
};

#include "TestCase1_def.hpp"

#endif  // TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_1_DECL_HPP
