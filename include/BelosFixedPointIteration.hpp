// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_FIXEDPOINT_ITERATION_HPP
#define BELOS_FIXEDPOINT_ITERATION_HPP

/*! \file BelosFixedPointIteration.hpp
    \brief Pure virtual base class which augments the basic interface for a fixed point linear solver iteration.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosIteration.hpp"

namespace Belos {

  //! @name FixedPointIteration Structures 
  //@{ 
  
  /** \brief Structure to contain pointers to FixedPointIteration state variables.
   *
   * This struct is utilized by FixedPointIteration::initialize() and FixedPointIteration::getState().
   */
  template <class ScalarType, class MV>
  struct FixedPointIterationState {

    /*! \brief The current residual. */
    Teuchos::RCP<const MV> R;

    /*! \brief The current preconditioned residual. */
    Teuchos::RCP<const MV> Z;
    
    FixedPointIterationState() : R(Teuchos::null), Z(Teuchos::null)
    {}
  };

template<class ScalarType, class MV, class OP>
class FixedPointIteration : virtual public Iteration<ScalarType,MV,OP> {

  public:

  //! @name State methods
  //@{ 
  /*! \brief Initialize the solver to an iterate, providing a complete state.
   *
   * The %FixedPointIteration contains a certain amount of state.
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
  virtual void initializeFixedPoint(FixedPointIterationState<ScalarType,MV>& newstate) = 0;

  /*! \brief Get the current state of the linear solver.
   *
   * The data is only valid if isInitialized() == \c true.
   *
   * \returns A FixedPointIterationState object containing const pointers to the current solver state.
   */
  virtual FixedPointIterationState<ScalarType,MV> getState() const = 0;
  //@}

};

} // end Belos namespace

#endif /* BELOS_FIXEDPOINT_ITERATION_HPP */
