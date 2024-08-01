//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_WrapperModelEvaluator_hpp
#define Tempus_WrapperModelEvaluator_hpp

#include "Tempus_config.hpp"
#include "Tempus_TimeDerivative.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace Tempus {

/// EVALUATION_TYPE indicates the evaluation to apply to the implicit ODE.
enum EVALUATION_TYPE {
  EVALUATE_RESIDUAL,      ///< Evaluate residual for the implicit ODE
  SOLVE_FOR_X,            ///< Solve for x and determine xDot from x.
  SOLVE_FOR_XDOT_CONST_X  ///< Solve for xDot keeping x constant (for ICs).
};

template <class Scalar>
class ImplicitODEParameters {
 public:
  /// Constructor
  ImplicitODEParameters()
    : timeDer_(Teuchos::null),
      timeStepSize_(Scalar(0.0)),
      alpha_(Scalar(0.0)),
      beta_(Scalar(0.0)),
      evaluationType_(SOLVE_FOR_X),
      stageNumber_(0)
  {
  }
  /// Constructor
  ImplicitODEParameters(Teuchos::RCP<TimeDerivative<Scalar> > timeDer,
                        Scalar timeStepSize, Scalar alpha, Scalar beta,
                        EVALUATION_TYPE evaluationType = SOLVE_FOR_X,
                        int stageNumber                = 0)
    : timeDer_(timeDer),
      timeStepSize_(timeStepSize),
      alpha_(alpha),
      beta_(beta),
      evaluationType_(evaluationType),
      stageNumber_(stageNumber)
  {
  }

  Teuchos::RCP<TimeDerivative<Scalar> > timeDer_;
  Scalar timeStepSize_;
  Scalar alpha_;
  Scalar beta_;
  EVALUATION_TYPE evaluationType_;
  int stageNumber_;
};

/** \brief A ModelEvaluator which wraps the application ModelEvaluator.
 *
 *  The WrapperModelEvaluator takes a state, \f$x\f$, computes time
 *  derivative(s), \f$\dot{x}\f$ and/or \f$\ddot{x}\f$, from the
 *  implicit stepper (StepperImplicit) and calls the application
 *  ModelEvaluator to determine its residual, \f$\mathcal{F}(x)\f$,
 *  which is suitable for the nonlinear solve.
 */
template <typename Scalar>
class WrapperModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar> {
 public:
  /// \name Vector Methods.
  //@{
  /// Get the x-solution space
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space()
      const = 0;

  /// Get the g space
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(
      int i) const = 0;

  /// Get the p space
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(
      int i) const = 0;
  //@}

  /// Set the underlying application ModelEvaluator
  virtual void setAppModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& me) = 0;

  /// Get the underlying application ModelEvaluator
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getAppModel()
      const = 0;

  /// Set parameters for application implicit ModelEvaluator solve.
  virtual void setForSolve(
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
      const Teuchos::RCP<ImplicitODEParameters<Scalar> >& p,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& y = Teuchos::null,
      const int index                                   = -1 /* index and y are for IMEX_RK_Partition */) = 0;
};

}  // namespace Tempus

#endif  // Tempus_WrapperModelEvaluator_hpp
