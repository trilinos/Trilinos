//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// ************************************************************************
//@HEADER

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
    ScalarType lambda;

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
  
  /** \brief LSQRIterationInitFailure is thrown when the LSQRIteration object is
   * unable to generate an initial iterate in the initialize() routine. 
   *
   * This std::exception is thrown from the initialize() method, which is
   * called by the user or from the iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this std::exception is thrown, 
   * isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
class LSQRIterationInitFailure : public BelosError {public:
      LSQRIterationInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief LSQRIterateFailure is thrown when the LSQRIteration object is unable to
   * compute the next iterate in the iterate() routine. 
   *
   * This std::exception is thrown from the iterate() method.
   *
   */
class LSQRIterateFailure : public BelosError {public:
      LSQRIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief LSQRIterationOrthoFailure is thrown when the LSQRIteration object is unable to
   * compute independent direction vectors in the iterate() routine. 
   *
   * This std::exception is thrown from the iterate() method.
   *
   */
class LSQRIterationOrthoFailure : public BelosError {public:
      LSQRIterationOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief LSQRIterationLAPACKFailure is thrown when a nonzero return value is passed back
   * from an LAPACK routine.
   *
   * This std::exception is thrown from the iterate() method.
   *
   */
class LSQRIterationLAPACKFailure : public BelosError {public:
      LSQRIterationLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}

} // end Belos namespace


#endif /* BELOS_LSQR_ITERATION_HPP */
