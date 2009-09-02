// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef BELOS_CG_ITERATION_HPP
#define BELOS_CG_ITERATION_HPP

/*! \file BelosCGIteration.hpp
    \brief Pure virtual base class which augments the basic interface for a conjugate gradient linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

namespace Belos {

  //! @name CGIteration Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to CGIteration state variables.
   *
   * This struct is utilized by CGIteration::initialize() and CGIteration::getState().
   */
  template <class ScalarType, class MV>
  struct CGIterationState {

    /*! \brief The current residual. */
    Teuchos::RCP<const MV> R;

    /*! \brief The current preconditioned residual. */
    Teuchos::RCP<const MV> Z;

    /*! \brief The current decent direction vector */
    Teuchos::RCP<const MV> P;

    /*! \brief The matrix A applied to current decent direction vector */
    Teuchos::RCP<const MV> AP;
    
    CGIterationState() : R(Teuchos::null), Z(Teuchos::null), 
		    P(Teuchos::null), AP(Teuchos::null)
    {}
  };

  //! @name CGIteration Exceptions
  //@{ 
  
  /** \brief CGIterationInitFailure is thrown when the CGIteration object is unable to
   * generate an initial iterate in the CGIteration::initialize() routine. 
   *
   * This std::exception is thrown from the CGIteration::initialize() method, which is
   * called by the user or from the CGIteration::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this std::exception is thrown, 
   * CGIteration::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the iteration.
   */
  class CGIterationInitFailure : public BelosError {public:
    CGIterationInitFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief CGIterateFailure is thrown when the CGIteration object is unable to
   * compute the next iterate in the CGIteration::iterate() routine. 
   *
   * This std::exception is thrown from the CGIteration::iterate() method.
   *
   */
  class CGIterateFailure : public BelosError {public:
    CGIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  /** \brief CGIterationOrthoFailure is thrown when the CGIteration object is unable to
   * compute independent direction vectors in the CGIteration::iterate() routine. 
   *
   * This std::exception is thrown from the CGIteration::iterate() method.
   *
   */
  class CGIterationOrthoFailure : public BelosError {public:
    CGIterationOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  /** \brief CGIterationLAPACKFailure is thrown when a nonzero return value is passed back
   * from an LAPACK routine.
   *
   * This std::exception is thrown from the CGIteration::iterate() method.
   *
   */
  class CGIterationLAPACKFailure : public BelosError {public:
    CGIterationLAPACKFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}


template<class ScalarType, class MV, class OP>
class CGIteration : virtual public Iteration<ScalarType,MV,OP> {

  public:

  //! @name State methods
  //@{ 
  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %CGIteration contains a certain amount of state, consisting of the current 
   * residual, preconditioned residual, and decent direction.
   *
   * initialize() gives the user the opportunity to manually set these,
   * although only the current unpreconditioned residual is required.
   *
   * \post 
   * <li>isInitialized() == \c true (see post-conditions of isInitialize())
   *
   * \note For any pointer in \c newstate which directly points to the multivectors in 
   * the solver, the data is not copied.
   */
  virtual void initializeCG(CGIterationState<ScalarType,MV> newstate) = 0;

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A CGIterationState object containing const pointers to the current solver state.
   */
  virtual CGIterationState<ScalarType,MV> getState() const = 0;
  //@}

};

} // end Belos namespace

#endif /* BELOS_CG_ITERATION_HPP */
