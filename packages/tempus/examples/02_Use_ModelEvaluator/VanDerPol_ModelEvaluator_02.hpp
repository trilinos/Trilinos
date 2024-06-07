//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_VANDERPOL_MODELEVALUATOR_02_HPP
#define TEMPUS_VANDERPOL_MODELEVALUATOR_02_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation


/** \brief ModelEvaluator implementation for the example \ref vanderpol.
  *
  *  This is a trimmed down version of the ModelEvaluator used in Tempus
  *  testing, Tempus_Test::VanDerPolModel, which exercises additional
  *  functionalities.  In this ModelEvaluator, we are just trying to
  *  demonstrate the mechanisms needed for the simple time integration,
  *  e.g., explicit and implicit ODEs for the \ref vanderpol.
  */
template<class Scalar>
class VanDerPol_ModelEvaluator_02
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
  public:

  /** Default Constructor sets the \ref vanderpol
   *  - parameter (\f$\epsilon=0.1\f$ )
   *  - the inital time (\f$t_0=0\f$ )
   *  - the initial conditions (\f$x_0(t_0=0) = 2\f$, \f$x_1(t_0=0) = 0\f$,
   *    \f$\dot{x}_0(t_0=0) = 0\f$, and \f$\dot{x}_1(t_0=0) = -2/\epsilon\f$ )
   *
   *  constructs
   *  - the solution vector space (a defaultSpmdVectorSpace of dimension 2)
   *  - the function evaluation vector space (a defaultSpmdVectorSpace of dimension 2)
   *  - a prototypical InArgs that just supports
   *    - the time (IN_ARGS_t)
   *    - the solution vector (IN_ARGS_x)
   *    - the time derivative of the solution vector (IN_ARGS_x_dot)
   *  - a prototypical OutArgs that just supports
   *    - the evaluation vector (OUT_ARG_f)
   *  - the nominal values (a.k.a. initial conditions for transient problems)
   *    - simply the prototypical InArgs with the initial condition vectors
   */
  VanDerPol_ModelEvaluator_02();

  /** \name Public functions overridden from ModelEvaluator. */
  //@{
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const
      { return x_space_; }
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const
      { return f_space_; }
    Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const
      { return nominalValues_; }
    Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const
      { return prototypicalInArgs_; }
  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const
      { return prototypicalOutArgs_; }
    /** Evaluate the model based on if the InArgs.get_x_dot() is null or
     *  not.  If it is null, evaluate the explicit ODE form of the model
     *  ( \f$ \dot{x} = f(x,t)\f$ ),
     *  \f{eqnarray*}{
     *    \dot{x}_0(t) & = & f_0 = x_1(t) \\
     *    \dot{x}_1(t) & = & f_1 = [(1-x_0^2)x_1-x_0]/\epsilon
     *  \f}
     *  otherwise, evaluate the implicit ODE form,
     *  \f$ \mathcal{F}(\dot{x},x,t) = 0 \f$,
     *  \f{eqnarray*}{
     *    \mathcal{F}_0 & = & \dot{x}_0(t) - x_1(t) \\
     *    \mathcal{F}_1 & = & \dot{x}_1(t) - [(1-x_0^2)x_1-x_0]/\epsilon
     *  \f}
     *  Both return their result through the OutArgs.get_f().
     */
    void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
      ) const;
  //@}

  /** Prototypical InArgs that just supports
   *  - the time (IN_ARGS_t)
   *  - the solution vector (IN_ARGS_x)
   *  - the time derivative of the solution vector (IN_ARGS_x_dot)
   */
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  prototypicalInArgs_;
  /// Prototypical OutArgs that just supports the evaluation vector (OUT_ARG_f)
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypicalOutArgs_;
  /** Nominal values (a.k.a. initial conditions for transient problems),
   *  and simply the prototypical InArgs with the initial condition vectors.
   */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  /// Solution vector space (a defaultSpmdVectorSpace of dimension 2)
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  /// Function evaluation vector space (a defaultSpmdVectorSpace of dimension 2)
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;

  int dim_;        ///< Number of state unknowns (2)
  Scalar t0_ic_;   ///< initial time = 0
  Scalar epsilon_; ///< This is a model parameter (\f$\epsilon=0.1\f$ )
  Scalar x0_ic_;   ///< initial condition for \f$x_0 = 2\f$
  Scalar x1_ic_;   ///< initial condition for \f$x_1 = 0\f$
};


#endif // TEMPUS_VANDERPOL_MODELEVALUATOR_02_HPP
