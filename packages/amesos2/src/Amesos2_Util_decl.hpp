// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

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
