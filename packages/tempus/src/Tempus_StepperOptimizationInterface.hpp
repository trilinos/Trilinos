//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Stepper_Optimization_Interface_hpp
#define Tempus_Stepper_Optimization_Interface_hpp

// Teuchos
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Tempus_config.hpp"

namespace Tempus {

/** \brief Stepper interface to support full-space optimization */
/*!
 * This is a potential interface to support transient full-space optimizations
 * methods such as those implemented by ROL through its DynamicConstraint
 * interface.  This interface is subject to major revision!
 *
 * Design consideration:  Take array of solution vectors as input, or
 * solution history, or array of states?
 *
 * The length of x is determined by the time step stencil, e.g., in an m-step
 * BDF method x would contain x_n,...,x_{n-m} at time step n in x[0],...,x[m]
 * with time values stored in t similarly.  p is the vector of design
 * parameters and param_index determines which model parameter this
 * corresponds to.
 */
template <class Scalar>
class StepperOptimizationInterface {
 public:
  StepperOptimizationInterface() {}

  virtual ~StepperOptimizationInterface() {}

  //! Return the number of solution vectors in the time step stencil
  virtual int stencilLength() const = 0;

  //! Compute time step residual
  virtual void computeStepResidual(
      Thyra::VectorBase<Scalar>& residual,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index) const = 0;

  //! Compute time step Jacobian
  /*!
   * deriv_index determines which component of x the derivative should be
   * computed with respect to.
   */
  virtual void computeStepJacobian(
      Thyra::LinearOpBase<Scalar>& jacobian,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index, const int deriv_index) const = 0;

  //! Compute time step derivative w.r.t. model parameters
  virtual void computeStepParamDeriv(
      Thyra::LinearOpBase<Scalar>& deriv,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index) const = 0;

  //! Compute time step Jacobian solver
  /*!
   * Derivative is always w.r.t. the most current solution vector
   */
  virtual void computeStepSolver(
      Thyra::LinearOpWithSolveBase<Scalar>& jacobian_solver,
      const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t, const Thyra::VectorBase<Scalar>& p,
      const int param_index) const = 0;
};

}  // namespace Tempus

#endif  // Tempus_Stepper_hpp
