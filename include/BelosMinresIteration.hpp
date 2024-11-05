// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_MINRES_ITERATION_HPP
#define BELOS_MINRES_ITERATION_HPP

/*! \file BelosMinresIteration.hpp
    \brief Pure virtual base class which augments the basic interface for a minimal residual linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

namespace Belos {

  //! @name MinresIteration Structures
  //@{ 
  
  /** \brief Structure to contain pointers to MinresIteration state variables.
   *
   * This struct is utilized by MinresIteration::initialize() and MinresIteration::getState().
   */
  template <class ScalarType, class MV>
  struct MinresIterationState {

    /*! \brief The current residual. */
    Teuchos::RCP<const MV> Y;

    /*! \brief Previous residual. */
    Teuchos::RCP<const MV> R1;

    /*! \brief Previous residual. */
    Teuchos::RCP<const MV> R2;

    /*! \brief The current direction vector */
    Teuchos::RCP<const MV> W;

    /*! \brief Previous direction vector */
    Teuchos::RCP<const MV> W1;

    /*! \brief Previous direction vector */
    Teuchos::RCP<const MV> W2;

    MinresIterationState() : Y(Teuchos::null),
                             R1(Teuchos::null),
                             R2(Teuchos::null),
                             W(Teuchos::null),
                             W1(Teuchos::null),
                             W2(Teuchos::null)
    {}
  };

  //! @name MinresIteration Exceptions
  //@{ 
  
  /** \brief MinresIterateFailure is thrown when the MinresIteration object is unable to
   * compute the next iterate in the MinresIteration::iterate() routine.
   *
   * This std::exception is thrown from the MinresIteration::iterate() method.
   *
   */
  class MinresIterateFailure : public BelosError {public:
    MinresIterateFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};
  
  //@}


template<class ScalarType, class MV, class OP>
class MinresIteration : virtual public Iteration<ScalarType,MV,OP> {

  public:

  //! @name State methods
  //@{ 
  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %MinresIteration contains a certain amount of state, consisting of the current
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
  virtual void initializeMinres(const MinresIterationState<ScalarType,MV> & newstate) = 0;

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A MinresIterationState object containing const pointers to the current solver state.
   */
  virtual MinresIterationState<ScalarType,MV> getState() const = 0;
  //@}

};

} // end Belos namespace

#endif /* BELOS_MINRES_ITERATION_HPP */
