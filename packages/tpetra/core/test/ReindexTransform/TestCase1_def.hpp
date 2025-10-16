// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_1_DEF_HPP
#define TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_1_DEF_HPP

#include "TestCase1_decl.hpp"
#include <MatrixMarket_Tpetra.hpp>

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t>::TestCase1(Teuchos::RCP<Teuchos::Comm<int> const> comm)
  : TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>(comm, 17  // globalNumRows
                                                          ,
                                                          17  // globalNumCols
                                                          ,
                                                          0  // rowIndexBase
                                                          ,
                                                          0  // colIndexBase
                                                          ,
                                                          5  // maxNnzPerRow
                                                          ,
                                                          43  // globalNumNnz
                                                          ,
                                                          4.940647893865e+00  // frobNorm
                                                          ,
                                                          {1.000000073665e+00, 1.000000073665e+00}  // lhsNorms2
                                                          ,
                                                          {1.000000075938e+00, 1.000000075938e+00}  // rhsNorms2
    ) {
  this->baseInstantiateLinearProblem("matA1.txt", "lhs1.txt", "rhs1.txt");

  this->baseCheckBeforeOrAfterTransform(this->m_linearProblem.get()  // Input
                                        ,
                                        "Case1_original"  // Input
  );

  this->m_mapsShallChange = false;
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t>::~TestCase1() {
  // Nothing to do
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
bool TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t>::checkTransformedProblem(Problem_t const* transformedProblem) const {
  return this->baseCheckTransformedProblem(transformedProblem);
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
bool TestCase1<Scalar_t, LocalId_t, GlobalId_t, Node_t>::checkAfterFwdRvs(Problem_t const* originalProblem) const {
  return this->baseCheckAfterFwdRvs(originalProblem);
}

#endif  // TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_1_DEF_HPP
