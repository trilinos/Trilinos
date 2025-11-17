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

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::TestCaseBase(Teuchos::RCP<Teuchos::Comm<int> const> comm, global_size_t const globalNumRows, global_size_t const globalNumCols, global_size_t const rowIndexBase, global_size_t const colIndexBase, global_size_t const maxNnzPerRow, global_size_t const globalNumNnz, Scalar_t const frobNorm, std::vector<Scalar_t> const &lhsNorms2, std::vector<Scalar_t> const &rhsNorms2)
  : m_comm(comm)
  , m_numRanks(m_comm->getSize())
  , m_myRank(m_comm->getRank())
  , m_globalNumRows(globalNumRows)
  , m_globalNumCols(globalNumCols)
  , m_rowIndexBase(rowIndexBase)
  , m_colIndexBase(colIndexBase)
  , m_maxNnzPerRow(maxNnzPerRow)
  , m_globalNumNnz(globalNumNnz)
  , m_frobNorm(frobNorm)
  , m_lhsNorms2(lhsNorms2)
  , m_rhsNorms2(rhsNorms2)
  , m_matrix(nullptr)
  , m_lhs(nullptr)
  , m_rhs(nullptr)
  , m_matrixRCP(Teuchos::null)
  , m_lhsRCP(Teuchos::null)
  , m_rhsRCP(Teuchos::null)
  , m_linearProblem(nullptr)
  , m_matrixMapsShallChange(false) {
  if (m_lhsNorms2.size() != m_rhsNorms2.size()) {
    std::stringstream msg;
    msg << "In TestCaseBase::constructor()"
        << ": vectors of norms should have equal sizes"
        << ", m_lhsNorms2.size() = " << m_lhsNorms2.size()
        << ", m_rhsNorms2.size() = " << m_rhsNorms2.size()
        << std::endl;
    throw std::runtime_error(msg.str());
  }
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::~TestCaseBase() {
  // Nothing to do
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
void TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseInstantiateLinearProblem(std::string const &matInputFileName, std::string const &lhsInputFileName, std::string const &rhsInputFileName) {
  // Create matrix
  Teuchos::RCP<Matrix_t> mat = Tpetra::MatrixMarket::Reader<Matrix_t>::readSparseFile(matInputFileName  // const std::string  & filename
                                                                                      ,
                                                                                      m_comm  // const trcp_tcomm_t & comm
                                                                                      ,
                                                                                      true  // const bool           callFillComplete
                                                                                      ,
                                                                                      false  // const bool           tolerant
                                                                                      ,
                                                                                      false  // const bool           debug
  );
  m_matrix                   = std::unique_ptr<Matrix_t>(new Matrix_t(*mat));
  m_matrixRCP                = Teuchos::rcp<Matrix_t>(m_matrix.get(), false);

  // Get a RCP to map
  Teuchos::RCP<Map_t const> mapRCP = Teuchos::rcp<Map_t const>(m_matrix->getRowMap().get(), false);

  // Create lhs
  Teuchos::RCP<MultiVector_t> lhs = Tpetra::MatrixMarket::Reader<MultiVector_t>::readDenseFile(lhsInputFileName  // const std::string            & filename
                                                                                               ,
                                                                                               m_comm  // const trcp_tcomm_t           & comm
                                                                                               ,
                                                                                               mapRCP  // Teuchos::RCP<const map_type> & map
                                                                                               ,
                                                                                               false  // const bool                     tolerant
                                                                                               ,
                                                                                               false  // const bool                     debug
                                                                                               ,
                                                                                               false  // const bool                     binary
  );
  m_lhs                           = std::unique_ptr<MultiVector_t>(new MultiVector_t(*lhs));
  m_lhsRCP                        = Teuchos::rcp<MultiVector_t>(m_lhs.get(), false);

  // Create rhs
  Teuchos::RCP<MultiVector_t> rhs = Tpetra::MatrixMarket::Reader<MultiVector_t>::readDenseFile(rhsInputFileName  // const std::string            & filename
                                                                                               ,
                                                                                               m_comm  // const trcp_tcomm_t           & comm
                                                                                               ,
                                                                                               mapRCP  // Teuchos::RCP<const map_type> & map
                                                                                               ,
                                                                                               false  // const bool                     tolerant
                                                                                               ,
                                                                                               false  // const bool                     debug
                                                                                               ,
                                                                                               false  // const bool                     binary
  );
  m_rhs                           = std::unique_ptr<MultiVector_t>(new MultiVector_t(*rhs));
  m_rhsRCP                        = Teuchos::rcp<MultiVector_t>(m_rhs.get(), false);

  // Create linear problem
  m_linearProblem = std::unique_ptr<Problem_t>(new Problem_t(m_matrixRCP, m_lhsRCP, m_rhsRCP));
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
typename TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::Problem_t *
TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::linearProblem() const {
  return m_linearProblem.get();
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
bool TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseCheckBeforeOrAfterTransform(Problem_t const *problem, std::string const &caseString, bool &solverMapTransformationWouldChangeMatrixMaps) const {
  bool result(true);
  Scalar_t relDiffThershold(1.e-12);

  // *****************************************************************
  // Check the matrix in 'problem'
  // *****************************************************************
  Matrix_t const *matrix = dynamic_cast<Matrix_t *>(problem->getMatrix().get());

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

  // Check globalNumRows
  if (matrix->getGlobalNumRows() != m_globalNumRows) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": globalNumRows should be " << m_globalNumRows
        << ", but it is " << matrix->getGlobalNumRows()
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Check globalNumCols
  if (matrix->getGlobalNumCols() != m_globalNumCols) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": globalNumCols should be " << m_globalNumCols
        << ", but it is " << matrix->getGlobalNumCols()
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Check maxNnzPerRow
  if (matrix->getGlobalMaxNumRowEntries() != m_maxNnzPerRow) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": maxNnzPerRow should be " << m_maxNnzPerRow
        << ", but it is " << matrix->getGlobalMaxNumRowEntries()
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Check globalNumNnz
  if (matrix->getGlobalNumEntries() != m_globalNumNnz) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": globalNumNnz should be " << m_globalNumNnz
        << ", but it is " << matrix->getGlobalNumEntries()
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  // Check Frobenius norm
  {
    Scalar_t frobNorm(matrix->getFrobeniusNorm());
    Scalar_t relDiff(std::abs(1. - frobNorm / m_frobNorm));
    if (relDiff > relDiffThershold) {
      std::stringstream msg;
      msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
          << ": frobNorm should be " << m_frobNorm
          << ", but it is " << frobNorm
          << ", causing relDiff = " << relDiff
          << " > threshold = " << relDiffThershold
          << std::endl;
      throw std::runtime_error(msg.str());
    }
  }

  // *****************************************************************
  // Check the lhs multivector
  // *****************************************************************
  if (problem->getLHS()->getNumVectors() != m_lhsNorms2.size()) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": problem->getLHS()->getNumVectors() = " << problem->getLHS()->getNumVectors()
        << ", m_lhsNorms2.size() = " << m_lhsNorms2.size()
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  if (problem->getLHS()->getGlobalLength() != m_globalNumCols) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": problem->getLHS()->getGlobalLength() = " << problem->getLHS()->getGlobalLength()
        << ", m_globalNumCols = " << m_globalNumCols
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  {
    size_t numVectors(problem->getLHS()->getNumVectors());
    std::vector<Scalar_t> allNorms(numVectors);
    Teuchos::ArrayView<Scalar_t> allNorms2(allNorms.data(), allNorms.size());
    problem->getLHS()->norm2(allNorms2);
    for (size_t v(0); v < numVectors; ++v) {
      if (m_lhsNorms2[v] == 0.) {
        if (allNorms2[v] > relDiffThershold) {
          std::stringstream msg;
          msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
              << "lhs, v = " << v
              << ": allNorms2[v] = " << allNorms2[v]
              << std::endl;
          throw std::runtime_error(msg.str());
        }
      } else {
        Scalar_t relDiff(std::abs(1. - allNorms2[v] / m_lhsNorms2[v]));
        if (relDiff > relDiffThershold) {
          std::stringstream msg;
          msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
              << "lhs, v = " << v
              << ": allNorms2[v] = " << allNorms2[v]
              << ", m_lhsNorms2[v] = " << m_lhsNorms2[v]
              << std::endl;
          throw std::runtime_error(msg.str());
        }
      }
    }
  }

  // *****************************************************************
  // Check the rhs multivector
  // *****************************************************************
  if (problem->getRHS()->getNumVectors() != m_rhsNorms2.size()) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": problem->getRHS()->getNumVectors() = " << problem->getRHS()->getNumVectors()
        << ", m_rhsNorms2.size() = " << m_rhsNorms2.size()
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  if (problem->getRHS()->getGlobalLength() != m_globalNumRows) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
        << ": problem->getRHS()->getGlobalLength() = " << problem->getRHS()->getGlobalLength()
        << ", m_globalNumRows = " << m_globalNumRows
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  {
    size_t numVectors(problem->getRHS()->getNumVectors());
    std::vector<Scalar_t> allNorms(numVectors);
    Teuchos::ArrayView<Scalar_t> allNorms2(allNorms.data(), allNorms.size());
    problem->getRHS()->norm2(allNorms2);
    for (size_t v(0); v < numVectors; ++v) {
      if (m_rhsNorms2[v] == 0.) {
        if (allNorms2[v] > relDiffThershold) {
          std::stringstream msg;
          msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
              << "rhs, v = " << v
              << ": allNorms2[v] = " << allNorms2[v]
              << std::endl;
          throw std::runtime_error(msg.str());
        }
      } else {
        Scalar_t relDiff(std::abs(1. - allNorms2[v] / m_rhsNorms2[v]));
        if (relDiff > relDiffThershold) {
          std::stringstream msg;
          msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
              << "rhs, v = " << v
              << ": allNorms2[v] = " << allNorms2[v]
              << ", m_rhsNorms2[v] = " << m_rhsNorms2[v]
              << std::endl;
          throw std::runtime_error(msg.str());
        }
      }
    }
  }

  {
    MultiVector_t tmpRhs(problem->getLHS()->getMap(), problem->getLHS()->getNumVectors(), true /* zeroOut */);
    m_matrix->apply(*problem->getLHS(), tmpRhs);

    // update(): this = beta*this + alpha*A
    tmpRhs.update(-1.  // alpha
                  ,
                  *m_rhs  // A
                  ,
                  1.  // beta
    );

    size_t numVectors(tmpRhs.getNumVectors());
    std::vector<Scalar_t> allNorms(numVectors);
    Teuchos::ArrayView<Scalar_t> allNorms2(allNorms.data(), allNorms.size());
    tmpRhs.norm2(allNorms2);
    for (size_t v(0); v < numVectors; ++v) {
      if (allNorms2[v] > relDiffThershold) {
        std::stringstream msg;
        msg << "In TestCaseBase::baseCheckBeforeOrAfterTransform()"
            << "tmpRhs - m_rhs, v = " << v
            << ": allNorms2[v] = " << allNorms2[v]
            << std::endl;
        throw std::runtime_error(msg.str());
      }
    }
  }

  return result;
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
bool TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseCheckTransformedProblem(Problem_t const *transformedProblem) const {
  bool result(true);

  bool solverMapTransformationWouldChangeMatrixMaps(false);
  result = baseCheckBeforeOrAfterTransform(transformedProblem  // Input
                                           ,
                                           "After transform"  // Input
                                           ,
                                           solverMapTransformationWouldChangeMatrixMaps  // Output
  );
  if (result == false) return result;

  // *****************************************************************
  // Check the matrix in 'transformedProblem'
  // *****************************************************************
  if (solverMapTransformationWouldChangeMatrixMaps == true) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckTransformedProblem()"
        << ", numRanks = " << this->m_numRanks
        << ": solverMapTransformationWouldChangeMatrixMaps should be 'false'"
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  if (m_matrixMapsShallChange == false) {
    // Both matrices (before and after the SolverMap transform) should have
    // equal maps, allowing us to call the CrsMatrix<>::add() method for
    // calculating (i) the difference between the matrices, and (ii) the
    // Frobenius norm of such difference.
    Matrix_t *tmpPtr1 = dynamic_cast<Matrix_t *>(m_linearProblem->getMatrix().get());
    Matrix_t *tmpPtr2 = dynamic_cast<Matrix_t *>(transformedProblem->getMatrix().get());

    result = this->baseCheckMatricesAreEqual(*tmpPtr1, *tmpPtr2);
    if (result == false) return result;
  }

  // *****************************************************************
  // Check the lhs multivector
  // *****************************************************************

  // No extra checks beyond the checks performed by 'baseCheckBeforeOrAfterTransform()'

  // *****************************************************************
  // Check the rhs multivector
  // *****************************************************************

  // No extra checks beyond the checks performed by 'baseCheckBeforeOrAfterTransform()'

  return result;
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
bool TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseCheckAfterFwdRvs(Problem_t const *originalProblem) const {
  bool result(true);

  bool solverMapTransformationWouldChangeMatrixMaps(false);
  result = baseCheckBeforeOrAfterTransform(m_linearProblem.get()  // Input
                                           ,
                                           "After FwdRvs"  // Input
                                           ,
                                           solverMapTransformationWouldChangeMatrixMaps  // Output
  );
  if (result == false) return result;

  // *****************************************************************
  // Check the matrix in 'm_linearProblem' (the problem that was
  // transformed).
  // *****************************************************************
  Matrix_t *tmpPtr1 = dynamic_cast<Matrix_t *>(originalProblem->getMatrix().get());
  Matrix_t *tmpPtr2 = dynamic_cast<Matrix_t *>(m_linearProblem->getMatrix().get());

  result = this->baseCheckMatricesAreEqual(*tmpPtr1, *tmpPtr2);
  if (result == false) return result;

  // *****************************************************************
  // Check the lhs multivector
  // *****************************************************************

  // No extra checks beyond the checks performed by 'baseCheckBeforeOrAfterTransform()'

  // *****************************************************************
  // Check the rhs multivector
  // *****************************************************************

  // No extra checks beyond the checks performed by 'baseCheckBeforeOrAfterTransform()'

  return result;
}

template <class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t>
bool TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>::baseCheckMatricesAreEqual(Matrix_t const &mat1, Matrix_t const &mat2) const {
  bool result(true);

  if ((*(mat1.getMap()) == *(mat2.getMap())) &&
      (*(mat1.getRowMap()) == *(mat2.getRowMap())) &&
      (*(mat1.getColMap()) == *(mat2.getColMap())) &&
      (*(mat1.getDomainMap()) == *(mat2.getDomainMap())) &&
      (*(mat1.getRangeMap()) == *(mat2.getRangeMap()))) {
    // Ok
  } else {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckMatricesAreEqual()"
        << ", numRanks = " << this->m_numRanks
        << ": matrices should have equal maps"
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  Matrix_t diff(mat1, Teuchos::Copy);
  Teuchos::RCP<Teuchos::ParameterList> params(Teuchos::null);
  Teuchos::RCP<Tpetra::RowMatrix<Scalar_t, LocalId_t, GlobalId_t, Node_t> > aux = diff.add(-1.0, mat2, 1.0, diff.getDomainMap(), diff.getRangeMap(), params);

  Matrix_t *tmpPtr3 = dynamic_cast<Matrix_t *>(aux.get());

  Scalar_t frobNorm_1(mat1.getFrobeniusNorm());
  Scalar_t frobNorm_2(mat2.getFrobeniusNorm());
  Scalar_t frobNorm_diff(tmpPtr3->getFrobeniusNorm());

  Scalar_t normThershold(1.e-12);
  if (frobNorm_diff > normThershold) {
    std::stringstream msg;
    msg << "In TestCaseBase::baseCheckMatricesAreEqual()"
        << ", numRanks = " << this->m_numRanks
        << ": mat2 should be equal to mat1"
        << ". Frobenius norms (mat1, mat2, diff)"
        << "= " << std::scientific << std::setprecision(12) << frobNorm_1
        << ", " << std::scientific << std::setprecision(12) << frobNorm_2
        << ", " << std::scientific << std::setprecision(12) << frobNorm_diff
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  return result;
}

#endif  // TPETRA_CORE_TEST_SOLVER_MAP_TRANSFORM_TEST_CASE_BASE_DEF_HPP
