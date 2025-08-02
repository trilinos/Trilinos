// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DECL_HPP
#define TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DECL_HPP

#include <Tpetra_LinearProblem.hpp>

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
class TestCaseBase
{
public:
  using Map_t         = Tpetra::Map          <          LocalId_t, GlobalId_t, Node_t>;
  using MultiVector_t = Tpetra::MultiVector  <Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Matrix_t      = Tpetra::CrsMatrix    <Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using Problem_t     = Tpetra::LinearProblem<Scalar_t, LocalId_t, GlobalId_t, Node_t>;
  using ColValPair_t  = std::pair<LocalId_t /*localColIndex*/, Scalar_t>;

  TestCaseBase() = delete;

  TestCaseBase( Teuchos::RCP<Teuchos::Comm<int> const> comm
              , GlobalId_t const globalNumRows
              , GlobalId_t const globalNumCols
              , GlobalId_t const rowIndexBase
              , GlobalId_t const colIndexBase
              , GlobalId_t const maxNnzPerRow
              , GlobalId_t const globalNumNnz
              , Scalar_t   const frobNorm
              );

  TestCaseBase( TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> const & source ) = delete;

  TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> & operator=( TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t> const & rhs ) = delete;

  virtual ~TestCaseBase();

  Problem_t * linearProblem() const;
      
  virtual bool checkTransformedProblem( Problem_t const * transformedProblem ) const = 0;
    
protected:
  void baseInstantiateLinearProblem( std::string const & inputFileName );

  bool internalBasicChecks_beforeOrAfterTransform( Problem_t   const * problem
                                                 , std::string const & caseString
                                                 , bool              & solverMapTransformationWouldChangeMatrixMaps
                                                 ) const;

  bool baseCheckTransformedProblem( Problem_t const * transformedProblem ) const;

  Teuchos::RCP< const Teuchos::Comm<int> > m_comm;
  int const m_numRanks;
  int const m_myRank;

  GlobalId_t const m_globalNumRows;
  GlobalId_t const m_globalNumCols;

  GlobalId_t const m_rowIndexBase;
  GlobalId_t const m_colIndexBase;

  GlobalId_t const m_maxNnzPerRow;
  GlobalId_t const m_globalNumNnz;
  Scalar_t   const m_frobNorm;

  std::unique_ptr< Map_t>          m_rowMap;
  std::unique_ptr< Map_t>          m_colMap;
  std::unique_ptr< Map_t>          m_domainMap;
  std::unique_ptr< Map_t>          m_rangeMap;
  std::unique_ptr< Matrix_t >      m_matrix;
  std::unique_ptr< MultiVector_t > m_lhs;
  std::unique_ptr< MultiVector_t > m_rhs;
  Teuchos::RCP< Matrix_t >         m_matrixRCP;
  Teuchos::RCP< MultiVector_t >    m_lhsRCP;
  Teuchos::RCP< MultiVector_t >    m_rhsRCP;
  std::unique_ptr< Problem_t >     m_linearProblem;
  bool                             m_mapsShallChange;
};

#include "TestCaseBase_def.hpp"

#endif // TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DECL_HPP
