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
  \file   Amesos2_Factory_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu Dec 31 10:32:08 2009

  \brief  Amesos2 Abstract Factory for producing solver interfaces. Declarations.
*/

#ifndef AMESOS2_FACTORY_DECL_HPP
#define AMESOS2_FACTORY_DECL_HPP

#include <Teuchos_RCP.hpp>

#include "Amesos2_config.h"

#include "Amesos2_SolverBase.hpp"
// include here klu2 headers

#ifdef HAVE_AMESOS2_KLU
#include "Amesos2_Klu.hpp"
#endif
#ifdef HAVE_AMESOS2_LAPACK
#include "Amesos2_Lapack.hpp"
#endif
#ifdef HAVE_AMESOS2_MUMPS
#include "Amesos2_Mumps.hpp"
#endif
#ifdef HAVE_AMESOS2_SCALAPACK
#include "Amesos2_Scalapack.hpp"
#endif
#ifdef HAVE_AMESOS2_UMFPACK
#include "Amesos2_Umfpack.hpp"
#endif
#ifdef HAVE_AMESOS2_SUPERLUDIST
#include "Amesos2_Superludist.hpp"
#endif
#ifdef HAVE_AMESOS2_SUPERLU
#include "Amesos2_Superlu.hpp"
#endif
#ifdef HAVE_AMESOS2_DSCPACK
#include "Amesos2_Dscpack.hpp"
#endif
#ifdef HAVE_AMESOS2_PARDISO
#include "Amesos2_Pardiso.hpp"
#endif
#ifdef HAVE_AMESOS2_TAUCS
#include "Amesos2_Taucs.hpp"
#endif
#ifdef HAVE_AMESOS2_PARAKLETE
#include "Amesos2_Paraklete.hpp"
#endif


using Teuchos::rcp;

namespace Amesos {


/**
 * \brief Abstract Factory for creating instances of Amesos2 Solver interfaces.
 *
 * \attention
 * Users should favour these static factory methods for creating Amesos2 solver
 * instances over explicitly instantiating their own.
 *
 * \note A solver's third-party library must be enabled in the Trilinos build,
 * and Amesos must be also told to enable it.  Put
 * <tt>Amesos2_ENABLE_<i>SOLVERNAME</i>:BOOL=ON</tt> in your Trilinos
 * configuration script to do this, where <i>SOLVERNAME</i> is the name of the
 * solver you would like to enable.
 *
 * \section usage Example Usage
 *
 * \code
 * typedef Tpetra::CrsMatrix<double,int> MAT;
 * typedef Tpetra::MultiVector<double,int> VEC;
 * // ... Create A of type RCP<MAT>, and X and B of type RCP<VEC> ...
 * RCP<SolverBase> solver = Factory<MAT,VEC>::create("Superlu", A, X, B);
 * \endcode
 *
 */
template <typename Matrix,
          typename Vector >
class Factory {
public:

  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * Default option which interfaces with the native KLU2 solver.
   *
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a KLU2 solver interface.
   */
  static Teuchos::RCP<SolverBase>
  create(Matrix* A, Vector* X, Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * Default option which interfaces with the native KLU2 solver.
   *
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a KLU2 solver interface.
   */
  static Teuchos::RCP<SolverBase>
  create(
    Teuchos::RCP<Matrix> A,
    Teuchos::RCP<Vector> X,
    Teuchos::RCP<Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName A \c C character string with the name of the
   * underlying third-party solver desired.
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   */
  static Teuchos::RCP<SolverBase>
  create(const char* solverName, Matrix* A, Vector* X, Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A pointer to the coefficient matrix
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   */
  static Teuchos::RCP<SolverBase>
  create(
    const char* solverName,
    const Teuchos::RCP<Matrix> A,
    const Teuchos::RCP<Vector> X,
    const Teuchos::RCP<Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   */
  static Teuchos::RCP<SolverBase>
  create(const std::string solverName, Matrix* A, Vector* X, Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   */
  static Teuchos::RCP<SolverBase>
  create(
    const std::string solverName,
    const Teuchos::RCP<Matrix> A,
    const Teuchos::RCP<Vector> X,
    const Teuchos::RCP<Vector> B);


  /**
   * \brief Queries the Factory for support of the named third-party library.
   *
   * \return \c true if the solver is supported.
   */
  static bool query(const char* solverName);


  /**
   * \brief Queries the Factory for support of the named third-party library.
   *
   * \return \c true if the solver is supported.
   */
  static bool query(const std::string solverName);


  // TODO: Here in Amesos.cpp there is a function defined
  // getValidParameters.  I wonder if it would be more appropriate
  // to define this function in the base Solver and concrete Solver
  // classes instead.

};				// end class Amesos2::Factory


} // end namespace Amesos

#endif	// AMESOS2_FACTORY_DECL_HPP
