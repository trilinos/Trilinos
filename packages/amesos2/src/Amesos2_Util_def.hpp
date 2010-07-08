/**
  \file   Amesos2_Utils_def.cpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed May 26 19:37:37 CDT 2010
  
  \brief  Utility functions for Amesos2
*/

#ifndef AMESOS2_UTIL_DEF_HPP
#define AMESOS2_UTIL_DEF_HPP

#include <Teuchos_ScalarTraits.hpp>

/** \brief Computes the true residual
 * 
 * Computes the true residual, B - A*X, and prints the results.
 *
 * \param A Matrix 
 * \param X Vector
 * \param B Vector
 * \param trans  \c true = use transpose for matrix-vector multiply
 * \param prefix string prefix to prints before rest of output
 *
 * TODO: Use Matrix and MultiVecAdapters instead of strait matrix and vector arguments
 */
template <typename Matrix,
          typename Vector>
void Amesos::Util::computeTrueResidual(
  const Teuchos::RCP<Matrix>& A,
  const Teuchos::RCP<Vector>& X,
  const Teuchos::RCP<Vector>& B,
  const Teuchos::ETransp trans=Teuchos::NO_TRANS,
  const std::string prefix="")
{
  // typename Teuchos::ScalarTraits<typename Matrix::scalar_type>::magnitudeType norm;
  // Tpetra::Vector<typename Matrix::scalar_type> Ax(B.getMap());
  // size_t numVectors = X.getNumVectors();

  // for( int i = 0; i < numVectors; ++i ){
  //   // TODO: Tpetra::Operator::apply only accepts
  //   // Tpetra::MultiVector objects, and indeed, only those that are
  //   // templated on the same types.  We may have to limit the
  //   // objects to Tpetra::RowMatrix and Tpetra::MultiVector
  //   // arguments.
  //   A.apply(*X.getVector(i), Ax, trans);
  //   Ax.update(1.0,*B.getVector(i),-1.0);
  //   norm = Ax.norm2();

  //   if (A.getComm().getRank() == 0){
  //     std::cout << prefix << " : vector "
  //               << i << ", ||b - Ax|| = "
  //               << norm << std::endl;
  //   }
  // }
}


/* We assume that Matrix and Vector are some instance of a
 * Amesos::MatrixAdapter or a Amesos::MultiVecAdapter, or at least implement
 * the required methods
 */
template< typename Matrix,
          typename Vector>
void Amesos::Util::computeVectorNorms(
  const Teuchos::RCP<Matrix> X,
  const Teuchos::RCP<Vector> B,
  std::string prefix="")
{
  typename Matrix::scalar_type normLHS, normRHS;
  size_t numVectors = X->getNumVectors();

  for (int i=0; i<numVectors; ++i){
    normLHS = X->getVector(i)->norm2();
    normRHS = B->getVector(i)->norm2();
    if (X->getMap()->getComm()->getRank() == 0){
      std::cout << prefix << " : vector "
                << ", ||x|| = " << normLHS
                << ", ||b|| = " << normRHS
                << std::endl;
    }
  }
}

template< typename Matrix>
void Amesos::Util::setMaxProcesses(
  const Teuchos::RCP<Matrix>& A,
  int& maxProcesses)
{
  int maxProcs = A->getComm()->getSize();

  switch(maxProcesses){
    case -3:
      maxProcesses = maxProcs; break;
    case -2:
      maxProcesses = (int) sqrt((double)maxProcs); break;
    case -1:			// We should do some testing on this
      // heuristic
      maxProcesses =
        1 + TEUCHOS_MAX(A.getGlobalNumRows() / 10000,
          A.getGlobalNumEntries() / 1000000);
      break;
  }

  if(maxProcesses <= 0) maxProcesses = 1;
  if(maxProcesses > maxProcs) maxProcesses = maxProcs;

  return;
}

/// Prints a line of 70 "-"s on std::cout.
void Amesos::Util::printLine() {
  std::cout << "-----------------------------------";
  std::cout << "-----------------------------------"
            << std::endl;
}
  

#endif	// #ifndef AMESOS2_UTIL_DEF_HPP
