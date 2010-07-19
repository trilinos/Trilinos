/**
  \file   Amesos2_Util_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu May 27 13:11:13 CDT 2010

  \brief  Utility functions for Amesos2
*/

#ifndef AMESOS2_UTIL_DECL_HPP
#define AMESOS2_UTIL_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_FancyOStream.hpp>

namespace Amesos {

namespace Util {


/**
 * \brief Computes the true residual
 *
 * Computes the true residual, \f$ B - A*X \f$ , and prints the results.
 *
 * \param A Matrix
 * \param X Vector
 * \param B Vector
 * \param trans  \c true = use transpose for matrix-vector multiply
 * \param prefix string prefix to prints before rest of output
 *
 */
template <typename Matrix,
          typename Vector>
void computeTrueResidual(
  const Teuchos::RCP<Matrix>& A,
  const Teuchos::RCP<Vector>& X,
  const Teuchos::RCP<Vector>& B,
  const Teuchos::ETransp trans=Teuchos::NO_TRANS,
  const std::string prefix="");


/* We assume that Matrix and Vector are some instance of a
 * Amesos::MatrixAdapter or a Amesos::MultiVecAdapter, or at least implement
 * the required methods
 */
template< typename Matrix,
          typename Vector>
void computeVectorNorms(
  const Teuchos::RCP<Matrix> X,
  const Teuchos::RCP<Vector> B,
  std::string prefix="");


/// Uses a heuristic to set the maximum number of processors
template <typename Matrix>
void setMaxProcesses(
  const Teuchos::RCP<Matrix>& A,
  int& maxProcesses);


/// Prints a line of 70 "-"s on std::cout.
void printLine( Teuchos::FancyOStream &out );


} // end namespace Util

} // end namespace Amesos

#endif	// #ifndef AMESOS2_UTIL_DECL_HPP
