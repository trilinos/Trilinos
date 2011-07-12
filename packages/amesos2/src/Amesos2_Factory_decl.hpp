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
 * \file   Amesos2_Factory_decl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun  7 16:25:16 2011
 * 
 * \brief  Contains declarations for Amesos::create and Amesos::query.
 *
 * Amesos2 defines the nonmember factory method \c Amesos::create for
 * creating instances of Amesos2 solvers.  The \c Amesos::query
 * function can be used to ask Amesos2 at runtime whether a particular
 * solver is supported.
 * 
 * \attention
 * Users should favor these factory methods for creating Amesos2 solver
 * instances over explicitly instantiating their own.
 *
 * \note A solver's third-party library must be enabled in the
 * Trilinos build, and Amesos2 must be also told to enable it.  Put
 * <tt>Amesos2_ENABLE_<i>SOLVERNAME</i>:BOOL=ON</tt> in your Trilinos
 * configuration script to do this, where <i>SOLVERNAME</i> is the
 * name of the solver you would like to enable.
 *
 * \section usage Example Usage
 *
 * \code
 * typedef Tpetra::CrsMatrix<double,int> MAT;
 * typedef Tpetra::MultiVector<double,int> VEC;
 * // ... Create A of type RCP<MAT>, and X and B of type RCP<VEC> ...
 * RCP<Amesos::Solver<MAT,VEC> > solver = Amesos::create<MAT,VEC>("SuperLU", A, X, B);
 * \endcode
 *
 * \ingroup amesos2_solver_framework
 */

#ifndef AMESOS2_FACTORY_DECL_HPP
#define AMESOS2_FACTORY_DECL_HPP

#include "Amesos2_config.h"

#ifdef HAVE_AMESOS2_KLU2
#include "Amesos2_Klu2.hpp"
#endif
// #ifdef HAVE_AMESOS2_MUMPS
// #include "Amesos2_Mumps.hpp"
// #endif
// #ifdef HAVE_AMESOS2_UMFPACK
// #include "Amesos2_Umfpack.hpp"
// #endif
#ifdef HAVE_AMESOS2_SUPERLUDIST	// Distributed-memory SuperLU
#include "Amesos2_Superludist.hpp"
#endif
#ifdef HAVE_AMESOS2_SUPERLUMT	// Multi-threaded SuperLU
#include "Amesos2_Superlumt.hpp"
#endif
#ifdef HAVE_AMESOS2_SUPERLU	// Sequential SuperLU
#include "Amesos2_Superlu.hpp"
#endif
// #ifdef HAVE_AMESOS2_DSCPACK
// #include "Amesos2_Dscpack.hpp"
// #endif
// #ifdef HAVE_AMESOS2_PARDISO
// #include "Amesos2_Pardiso.hpp"
// #endif
// #ifdef HAVE_AMESOS2_TAUCS
// #include "Amesos2_Taucs.hpp"
// #endif


namespace Amesos {

  template <class,class> class Solver;

  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * This is the default option which interfaces with the native KLU2 solver.
   *
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A \c C pointer to a KLU2 solver interface.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
  create(Matrix* A, Vector* X, Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * This is the default option which interfaces with the native KLU2 solver.
   *
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   * 
   * \return A <tt>Teuchos::RCP</tt> to a KLU2 solver interface.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(Teuchos::RCP<Matrix> A,
	 Teuchos::RCP<Vector> X,
	 Teuchos::RCP<Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName A \c C character string with the name of the
   *                        underlying third-party solver desired.
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A \c C pointer to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
  create(const char* solverName, Matrix* A, Vector* X, Vector* B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const char* solverName,
	 const Teuchos::RCP<Matrix> A,
	 const Teuchos::RCP<Vector> X,
	 const Teuchos::RCP<Vector> B);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A \c C pointer to the coefficient matrix
   * \param [in] X \c C pointer to LHS solution vector
   * \param [in] B \c C pointer to RHS vector
   *
   * \return A \c C pointer to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
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
   * \return A <tt>Teuchos::RCP</tt> to an Amesos2 solver interface.
   *
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string solverName,
	 const Teuchos::RCP<Matrix> A,
	 const Teuchos::RCP<Vector> X,
	 const Teuchos::RCP<Vector> B);

  
  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A.
   *
   * Suitable for cases where numeric factorization must be performed
   * before the X and B vectors are known.  Before a solve, the \c
   * setX() and \c setB() functions should be used to set X and B, or
   * the overloaded solve(X,B) method should be used.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A \c C pointer to the coefficient matrix
   *
   * \return A \c C pointer to an Amesos2 solver interface.
   * 
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Solver<Matrix,Vector>*
  create(const std::string solverName, Matrix* A);


  /**
   * \brief Creates an Amesos2 Solver interface with Matrix A.
   *
   * Suitable for cases where numeric factorization must be performed
   * before the X and B vectors are known.  Before a solve, the \c
   * setX() and \c setB() functions should be used to set X and B, or
   * the overloaded solve(X,B) method should be used.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to the coefficient matrix
   *
   * \return A <tt>Teuchos::RCP</tt> to an Amesos2 solver interface.
   * 
   * \throw std::invalid_argument The third-party solver named by \c
   * solverName is not supported.
   *
   * \relatesalso Amesos::Solver
   */
  template <class Matrix,
	    class Vector >
  Teuchos::RCP<Solver<Matrix,Vector> >
  create(const std::string solverName,
	 const Teuchos::RCP<Matrix> A);


  /**
   * \brief Queries the Factory for support of the named third-party library.
   *
   * \return \c true if the solver is supported.
   *
   * \relatesalso Amesos::Solver
   */
  bool query(const char* solverName);


  /**
   * \brief Queries the Factory for support of the named third-party library.
   *
   * \return \c true if the solver is supported.
   *
   * \relatesalso Amesos::Solver
   */
  bool query(const std::string solverName);

}

#endif	// AMESOS2_FACTORY_DECL_HPP
