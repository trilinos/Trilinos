// @HEADER
// ***********************************************************************
//
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/**
  \file   Amesos2_Factory_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu Dec 31 10:32:08 2009
  
  \brief Amesos2 Abstract Factory for producing solver interfaces.
*/

#ifndef AMESOS2_FACTORY_DECL_HPP
#define AMESOS2_FACTORY_DECL_HPP

#include "Amesos2_config.h"

#include "Amesos2_SolverBase.hpp"

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


#include <Teuchos_RCP.hpp>
using Teuchos::rcp;

namespace Amesos {


template <typename Matrix,
          typename Vector >
class Factory {
public:
  /* \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
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
   */
  static Teuchos::RCP<SolverBase>
  create(const char* solverName, Matrix* A, Vector* X, Vector* B);

  /* \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver 
   * \param [in] A pointer to a matrix of coefficients
   * \param [in] X pointer to LHS solution vector
   * \param [in] B pointer to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   */
  static Teuchos::RCP<SolverBase>
  create(
    const char* solverName,
    const Teuchos::RCP<Matrix> A,
    const Teuchos::RCP<Vector> X,
    const Teuchos::RCP<Vector> B);

  /* \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver 
   * \param [in] A <tt>Teuchos::RCP</tt> to a matrix of coefficients
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   */
  static Teuchos::RCP<SolverBase>
  create(const std::string solverName, Matrix* A, Vector* X, Vector* B);
  
  /* \brief Creates an Amesos2 Solver interface with Matrix A, LHS vector X,
   * and RHS vector B.
   *
   * \param [in] solverName The name of the desired third-party solver
   * \param [in] A <tt>Teuchos::RCP</tt> to a matrix of coefficients
   * \param [in] X <tt>Teuchos::RCP</tt> to LHS solution vector
   * \param [in] B <tt>Teuchos::RCP</tt> to RHS vector
   *
   * \return A <tt>Teuchos::RCP</tt> to a Amesos2 solver interface.
   */
  static Teuchos::RCP<SolverBase>
  create(
    const std::string solverName,
    const Teuchos::RCP<Matrix> A,
    const Teuchos::RCP<Vector> X,
    const Teuchos::RCP<Vector> B);

  /* \brief Queries the Factory for support of the named third-party library.
   */
  bool query(const char* solverName);

  /* \brief Queries the Factory for support of the named third-party library.
   */
  bool query(const std::string solverName);


  // TODO: Here in Amesos.cpp there is a function defined
  // getValidParameters.  I wonder if it would be more appropriate
  // to define this function in the base Solver and concrete Solver
  // classes instead.
    
};				// end class Amesos2::Factory

  
} // end namespace Amesos

#endif	// AMESOS2_FACTORY_DECL_HPP
