// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_CG_STOCHASTIC_ITERATION_HPP
#define BELOS_CG_STOCHASTIC_ITERATION_HPP

/*! \file BelosStochasticCGIteration.hpp
    \brief Pure virtual base class which augments the basic interface for a stochastic conjugate gradient linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"
#include "BelosCGIteration.hpp"

namespace Belos {

  //! @name CGIteration Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to CGIteration state variables.
   *
   * This struct is utilized by StochasticCGIteration::initialize() and StochasticCGIteration::getState().
   */
  template <class ScalarType, class MV>
  struct StochasticCGIterationState {

    /*! \brief The current residual. */
    Teuchos::RCP<const MV> R;

    /*! \brief The current preconditioned residual. */
    Teuchos::RCP<const MV> Z;

    /*! \brief The current decent direction vector */
    Teuchos::RCP<const MV> P;

    /*! \brief The matrix A applied to current decent direction vector */
    Teuchos::RCP<const MV> AP;
    
    /*! \brief The current stochastic recurrence vector */
    Teuchos::RCP<const MV> Y;

    StochasticCGIterationState() : R(Teuchos::null), Z(Teuchos::null), 
			 P(Teuchos::null), AP(Teuchos::null), Y(Teuchos::null)
    {}
  };

template<class ScalarType, class MV, class OP>
class StochasticCGIteration : virtual public Iteration<ScalarType,MV,OP> {

  public:

  //! @name State methods
  //@{ 
  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %StochasticCGIteration contains a certain amount of state, consisting of the current 
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
  virtual void initializeCG(StochasticCGIterationState<ScalarType,MV>& newstate) = 0;

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A StochasticCGIterationState object containing const pointers to the current solver state.
   */
  virtual StochasticCGIterationState<ScalarType,MV> getState() const = 0;
  //@}

};

} // end Belos namespace

#endif /* BELOS_STOCHASTIC_CG_ITERATION_HPP */
