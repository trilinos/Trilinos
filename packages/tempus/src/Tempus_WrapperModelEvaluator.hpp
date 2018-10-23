// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_WrapperModelEvaluator_hpp
#define Tempus_WrapperModelEvaluator_hpp

#include "Tempus_TimeDerivative.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace Tempus {

/// EVALUATION_TYPE indicates the evaluation to apply to the implicit ODE.
enum EVALUATION_TYPE {
  EVALUATE_RESIDUAL,      ///< Evaluate residual for the implicit ODE
  SOLVE_FOR_X,            ///< Solve for x and determine xDot from x.
  SOLVE_FOR_XDOT_CONST_X  ///< Solve for xDot keeping x constant (for ICs).
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
class WrapperModelEvaluator : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /// \name Vector Methods.
  //@{
    /// Get the x-solution space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_x_space() const = 0;

    /// Get the g space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_g_space(int i) const = 0;

    /// Get the p space
    virtual Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
      get_p_space(int i) const = 0;
  //@}

  /// Set the underlying application ModelEvaluator
  virtual void setAppModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me) = 0;

  /// Get the underlying application ModelEvaluator
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    getAppModel() const = 0;

  /// Set InArgs the wrapper ModelEvalutor.
  virtual void setInArgs(Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs) = 0;

  /// Get InArgs the wrapper ModelEvalutor.
  virtual Thyra::ModelEvaluatorBase::InArgs<Scalar> getInArgs() = 0;

  /// Set OutArgs the wrapper ModelEvalutor.
  virtual void setOutArgs(Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs)=0;

  /// Get OutArgs the wrapper ModelEvalutor.
  virtual Thyra::ModelEvaluatorBase::OutArgs<Scalar> getOutArgs() = 0;

  /// Set parameters for application implicit ModelEvaluator solve.
  virtual void setForSolve(Teuchos::RCP<TimeDerivative<Scalar> > td,
    Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs,
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs,
    EVALUATION_TYPE evaluationType = SOLVE_FOR_X) = 0;
};

} // namespace Tempus

#endif // Tempus_WrapperModelEvaluator_hpp
