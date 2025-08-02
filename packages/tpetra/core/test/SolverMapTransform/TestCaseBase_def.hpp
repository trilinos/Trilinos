// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DEF_HPP
#define TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DEF_HPP

#include "TestCaseBase_decl.hpp"
#include <MatrixMarket_Tpetra.hpp>

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::TestCaseBase( Teuchos::RCP<Teuchos::Comm<int> const> comm
                                                                   , global_size_t const globalNumRows
                                                                   , global_size_t const globalNumCols
                                                                   , global_size_t const rowIndexBase
                                                                   , global_size_t const colIndexBase
                                                                   , global_size_t const maxNnzPerRow
                                                                   , global_size_t const globalNumNnz
                                                                   , Scalar_t      const frobNorm
                                                                   )
  : m_comm           ( comm )
  , m_numRanks       ( m_comm->getSize() )
  , m_myRank         ( m_comm->getRank() )
  , m_globalNumRows  ( globalNumRows )
  , m_globalNumCols  ( globalNumCols )
  , m_rowIndexBase   ( rowIndexBase )
  , m_colIndexBase   ( colIndexBase )
  , m_maxNnzPerRow   ( maxNnzPerRow )
  , m_globalNumNnz   ( globalNumNnz )
  , m_frobNorm       ( frobNorm )
  , m_rowMap         ( nullptr )
  , m_colMap         ( nullptr )
  , m_domainMap      ( nullptr )
  , m_rangeMap       ( nullptr )
  , m_matrix         ( nullptr )
  , m_lhs            ( nullptr )
  , m_rhs            ( nullptr )
  , m_matrixRCP      ( Teuchos::null )
  , m_lhsRCP         ( Teuchos::null )
  , m_rhsRCP         ( Teuchos::null )
  , m_linearProblem  ( nullptr )
  , m_mapsShallChange( false )
{
  // Nothing to do
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::~TestCaseBase()
{
  // Nothing to do
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
void
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseInstantiateLinearProblem( std::string const & inputFileName )
{
  Teuchos::RCP<Matrix_t> mat = Tpetra::MatrixMarket::Reader<Matrix_t>::readSparseFile( inputFileName // const std::string  & filename
                                                                                     , m_comm        // const trcp_tcomm_t & comm
                                                                                     , true          // const bool           callFillComplete=true
                                                                                     , false         // const bool           tolerant=false
                                                                                     , false         // const bool           debug=false
                                                                                     );
  m_matrix = std::unique_ptr<Matrix_t>( new Matrix_t(*mat) );
  
  m_matrixRCP = Teuchos::rcp<Matrix_t     >(m_matrix.get(), false);
  m_lhsRCP    = Teuchos::rcp<MultiVector_t>(m_lhs.get(), false);
  m_rhsRCP    = Teuchos::rcp<MultiVector_t>(m_rhs.get(), false);
  m_linearProblem = std::unique_ptr<Problem_t>(new Problem_t(m_matrixRCP, m_lhsRCP, m_rhsRCP));
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::Problem_t *
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::linearProblem() const
{
  return m_linearProblem.get();
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
bool
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::internalBasicChecks_beforeOrAfterTransform
  ( Problem_t   const * problem
  , std::string const & caseString
  , bool              & solverMapTransformationWouldChangeMatrixMaps
  ) const
{
  bool result(true);
  
  Matrix_t const * matrix = dynamic_cast<Matrix_t*>(problem->getMatrix().get());

  size_t domainMap_localSize = matrix->getDomainMap()->getLocalNumElements();
  size_t localNumDifferences(0);
  for (size_t i(0); i < domainMap_localSize; ++i) {
    if (matrix->getDomainMap()->getGlobalElement(i) != matrix->getColMap()->getGlobalElement(i)) {
      localNumDifferences += 1;
      break;
    }
  }
  size_t globalNumDifferences(0);
  Teuchos::reduceAll(*m_comm, Teuchos::REDUCE_SUM, 1, &localNumDifferences, &globalNumDifferences);
  
  solverMapTransformationWouldChangeMatrixMaps = (globalNumDifferences != 0);
  if (false) { // if (solverMapTransformationWouldChangeMatrixMaps) {
    if (this->m_myRank == 0) {
      std::cout.flush();
      std::cout << caseString
                << ", numRanks = " << this->m_numRanks
                << ": solverMapTransformationWouldChangeMatrixMaps = " << solverMapTransformationWouldChangeMatrixMaps
                << std::endl;
      std::cout.flush();
    }
    this->m_comm->barrier();

    std::cout.flush();
    m_comm->barrier();
    for (int p(0); p < m_numRanks; ++p) {
      m_comm->barrier();
      if (p != m_myRank) continue;

      std::cout.flush();

      std::cout << caseString
                << ": rowMap()->globalNumElems = " << matrix->getRowMap()->getGlobalNumElements()
                << ", localNumElems = " << matrix->getRowMap()->getLocalNumElements()
                << "; globalIndices =";
      for (size_t i(0); i < matrix->getRowMap()->getLocalNumElements(); ++i) {
        std::cout << " " << matrix->getRowMap()->getGlobalElement(i);
      }
      std::cout << std::endl;

      std::cout << caseString
                << ": colMap()->globalNumElems = " << matrix->getColMap()->getGlobalNumElements()
                << ", localNumElems = " << matrix->getColMap()->getLocalNumElements()
                << "; globalIndices =";
      for (size_t i(0); i < matrix->getColMap()->getLocalNumElements(); ++i) {
        std::cout << " " << matrix->getColMap()->getGlobalElement(i);
      }
      std::cout << std::endl;

      std::cout << caseString
                << ": domainMap()->globalNumElems = " << matrix->getDomainMap()->getGlobalNumElements()
                << ", localNumElems = " << matrix->getDomainMap()->getLocalNumElements()
                << "; globalIndices =";
      for (size_t i(0); i < matrix->getDomainMap()->getLocalNumElements(); ++i) {
        std::cout << " " << matrix->getDomainMap()->getGlobalElement(i);
      }
      std::cout << std::endl;

      std::cout << caseString
                << ": rangeMap()->globalNumElems = " << matrix->getRangeMap()->getGlobalNumElements()
                << ", localNumElems = " << matrix->getRangeMap()->getLocalNumElements()
                << "; globalIndices =";
      for (size_t i(0); i < matrix->getRangeMap()->getLocalNumElements(); ++i) {
        std::cout << " " << matrix->getRangeMap()->getGlobalElement(i);
      }
      std::cout << std::endl;

      std::cout.flush();
    }
    m_comm->barrier();
  } // if (printDetailedInfo)

  // Check globalNumRows
  if (matrix->getGlobalNumRows() != m_globalNumRows) {
    std::stringstream msg;
    msg << "In TestCaseBase::internalBasicChecks_beforeOrAfterTransform()"
        << ": globalNumRows should be " << m_globalNumRows
        << ", but it is " << matrix->getGlobalNumRows()
        << std::endl;
    throw std::runtime_error( msg.str() );
  }

  // Check globalNumCols
  if (matrix->getGlobalNumCols() != m_globalNumCols) {
    std::stringstream msg;
    msg << "In TestCaseBase::internalBasicChecks_beforeOrAfterTransform()"
        << ": globalNumCols should be " << m_globalNumCols
        << ", but it is " << matrix->getGlobalNumCols()
        << std::endl;
    throw std::runtime_error( msg.str() );
  }

  // Check maxNnzPerRow
  if (matrix->getGlobalMaxNumRowEntries() != m_maxNnzPerRow) {
    std::stringstream msg;
    msg << "In TestCaseBase::internalBasicChecks_beforeOrAfterTransform()"
        << ": maxNnzPerRow should be " << m_maxNnzPerRow
        << ", but it is " << matrix->getGlobalMaxNumRowEntries()
        << std::endl;
    throw std::runtime_error( msg.str() );
  }

  // Check globalNumNnz
  if (matrix->getGlobalNumEntries() != m_globalNumNnz) {
    std::stringstream msg;
    msg << "In TestCaseBase::internalBasicChecks_beforeOrAfterTransform()"
        << ": globalNumNnz should be " << m_globalNumNnz
        << ", but it is " << matrix->getGlobalNumEntries()
        << std::endl;
    throw std::runtime_error( msg.str() );
  }

  // Check Frobenius norm
  Scalar_t relDiffThershold(1.e-12);
  Scalar_t frobNorm( matrix->getFrobeniusNorm() );
  Scalar_t relDiff( abs(1. - frobNorm/m_frobNorm) );
  if (relDiff > relDiffThershold) {
    std::stringstream msg;
    msg << "In TestCaseBase::internalBasicChecks_beforeOrAfterTransform()"
        << ": frobNorm should be " << m_frobNorm
        << ", but it is "          << frobNorm
        << ", causing relDiff = "  << relDiff
        << " > threshold = "       << relDiffThershold
        << std::endl;
    throw std::runtime_error( msg.str() );
  }

  return result;
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
bool
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseCheckTransformedProblem( Problem_t const * transformedProblem ) const
{
  bool result(true);

  bool solverMapTransformationWouldChangeMatrixMaps(false);
  result = internalBasicChecks_beforeOrAfterTransform( transformedProblem
                                                     , "After transform"
                                                     , solverMapTransformationWouldChangeMatrixMaps
                                                     );
  if (result == false) return result;

  if (this->m_myRank == 0) {
    std::cout.flush();
    std::cout << "After transform"
              << ", numRanks = " << this->m_numRanks
              << ", m_mapsShallChange = " << m_mapsShallChange
              << ": solverMapTransformationWouldChangeMatrixMaps = " << solverMapTransformationWouldChangeMatrixMaps
              << std::endl;
    std::cout.flush();
  }
  this->m_comm->barrier();

  if (solverMapTransformationWouldChangeMatrixMaps == true) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckTransformedProblem()"
        << ", numRanks = " << this->m_numRanks
        << ": solverMapTransformationWouldChangeMatrixMaps should be 'false'"
        << std::endl;
    throw std::runtime_error( msg.str() );
  }

  if (m_mapsShallChange == false) {
    Matrix_t * tmpPtr1 = dynamic_cast<Matrix_t *>( m_linearProblem->getMatrix().get() );
    Scalar_t frobNorm_1( tmpPtr1->getFrobeniusNorm() );

    Matrix_t * tmpPtr2 = dynamic_cast<Matrix_t *>( transformedProblem->getMatrix().get() );
    Scalar_t frobNorm_2( tmpPtr2->getFrobeniusNorm() );

    Matrix_t diff( *tmpPtr1, Teuchos::Copy );

    Teuchos::RCP<Teuchos::ParameterList> params(Teuchos::null);
    Teuchos::RCP< Tpetra::RowMatrix<Scalar_t, LocalId_t, GlobalId_t, Node_t> > aux = diff.add( -1.0
                                                                                             , *tmpPtr2
                                                                                             , 1.0
                                                                                             , diff.getDomainMap()
                                                                                             , diff.getRangeMap()
                                                                                             , params
                                                                                             );

    Matrix_t * tmpPtr3 = dynamic_cast<Matrix_t *>( aux.get() );
    Scalar_t frobNorm_diff( tmpPtr3->getFrobeniusNorm() );

    Scalar_t normThershold(1.e-12);
    if (frobNorm_diff > normThershold) {
      std::stringstream msg;
      msg << "In TestCaseBase::baseCheckTransformedProblem()"
          << ", numRanks = " << this->m_numRanks
          << ": because no map differences were detected in the original matrix, the matrix after transform should be equal to the original one"
          << ". Frobenius norms (orig, xform, diff)"
          << "= " << std::scientific << std::setprecision(12) << frobNorm_1
          << ", " << std::scientific << std::setprecision(12) << frobNorm_2
          << ", " << std::scientific << std::setprecision(12) << frobNorm_diff
          << std::endl;
      throw std::runtime_error( msg.str() );
    }
  }

  return result;
}
  

#endif // TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DEF_HPP
