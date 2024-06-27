// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_LSQR_ITERATION_HPP
#define BELOS_LSQR_ITERATION_HPP

/*! \file BelosLSQRIteration.hpp
    \brief IterationState contains the data that defines the state of
           the LSQR solver at any given time.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

namespace Belos {

  //! @name LSQRIteration Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to LSQRIteration state variables, ...
   *
   * This struct is utilized by initialize() and getState().
   * augment the basic interface for a Gmres linear solver iteration.
   */
  template <class ScalarType, class MV>
  struct LSQRIterationState {

    /*! \brief Bidiagonalization vector. */
    Teuchos::RCP<const MV> U;

    /*! \brief Bidiagonalization vector. */
    Teuchos::RCP<const MV> V;

    /*! \brief The search direction vector. */
    Teuchos::RCP<const MV> W;

    /*! \brief The damping value. */
    typename Teuchos::ScalarTraits<ScalarType>::magnitudeType lambda;

    /*! \brief The current residual norm. */
    ScalarType resid_norm;

    /*! \brief An approximation to the Frobenius norm of A. */
    ScalarType frob_mat_norm;

    /*! \brief An approximation to the condition number of A. */
    ScalarType mat_cond_num;

    /*! \brief An estimate of the norm of A^T*resid. */
    ScalarType mat_resid_norm;

    /*! \brief An estimate of the norm of the solution. */
    ScalarType sol_norm;

    /*! \brief The norm of the RHS vector b. */
    ScalarType bnorm;
    
    LSQRIterationState() : U(Teuchos::null), V(Teuchos::null), 
			   W(Teuchos::null), lambda(0.0), 
			   resid_norm(0.0), frob_mat_norm(0.0),
			   mat_cond_num(0.0), mat_resid_norm(0.0),
			   sol_norm(0.0), bnorm(0.0)
    {}
  };

  //! @name LSQRIteration Exceptions
  //@{ 
  
  /** \brief LSQRIterateFailure is thrown when the LSQRIteration object is unable to
   * compute the next iterate in the iterate() routine. 
   *
   * This std::exception is thrown from the iterate() method.
   *
   */
class LSQRIterateFailure : public BelosError {public:
      LSQRIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}

} // end Belos namespace


#endif /* BELOS_LSQR_ITERATION_HPP */
