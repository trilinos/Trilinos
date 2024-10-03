//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_HARMONIC_OSCILLATOR_MODEL_DECL_HPP
#define TEMPUS_TEST_HARMONIC_OSCILLATOR_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief
 * Consider the ODE:
 *   \f[
 *   m\ddot{x} + c\dot{x} + kx=f
 *   \f]
 * where \f$k \geq 0\f$ is a constant, \f$c\f$ is a constant damping parameter,
 * \f$f\f$ is a constant forcing parameter, and \f$m>0\f$ is a constant mass
 * parameter, with initial conditions are: \f{eqnarray*}{
 *     x(0) & = & 0\\
 *     \dot{x}(0) & = & 1
 *   \f}
 * It is straight-forward to show that the exact solution to this ODE is:
 *   \f{eqnarray*}{
 *     x(t) & = & t(1+0.5\tilde{f}t),  \hspace{3.6cm} if  \hspace{0.2cm} k = c =
 * 0 \\
 *          & = & \frac{(\tilde{c}-\tilde{f})}{\tilde{c}^2}(1-e^{-\tilde{c}t}) +
 * \frac{f}{c}t, \hspace{1.9cm} if \hspace{0.2cm} k = 0, c\neq 0 \\ & = &
 * \frac{1}{\sqrt{\tilde{k}}}\sin(\sqrt{\tilde{k}}t) +
 * \frac{f}{k}\left(1-\cos(\sqrt{\tilde{k}}t) \right), \hspace{0.2cm} if
 * \hspace{0.2cm} k > 0, c = 0 \f} where \f$\tilde{c}\equiv c/m\f$,
 * \f$\tilde{k}\equiv k/m\f$ and \f$\tilde{f}\equiv f/m\f$.  While it is
 * possible to derive the solution to this ODE for the case when \f$k > 0\f$ and
 * \f$c \neq 0\f$, we do not consider that case here.  When \f$c = k = 0\f$,
 * \f$m=1\f$, and \f$f=-1\f$, our ODE simplies to a canonical differential
 * equation model of a ball thrown up in the air, with a parabolic trajectory
 * solution, namely \f[ x(t) = t(1-0.5t) \f] where \f$t\in [0,2]\f$.  An
 * EpetraExt version of this simplified version of the test is implemented in
 * Piro::MockModelEval_B (see Trilinos/packages/piro/test), where it is used to
 * test the Piro (EpetraExt) Newmark-Beta scheme (see input_Solver_NB.xml input
 * file).  When \f$c = f = 0\f$ and  \f$m=k = 1\f$, this test is equivalent to
 * the SinCos model..
 */

template <class Scalar>
class HarmonicOscillatorModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase {
 public:
  // Constructor
  HarmonicOscillatorModel(
      Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null,
      const bool use_accel_IC                    = false);

  // Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  //@}

  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const &paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

 private:
  void setupInOutArgs_() const;

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar) const;
  //@}

 private:
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_vec_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_vec_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_dot_vec_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > p_init_;
  int vecLength_;     // Number of state unknowns (1)
  int numResponses_;  // Number of responses (1)
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  mutable bool isInitialized_;
  double c_;  // damping coefficient
  double f_;  // forcing coefficient
  double k_;  // coefficient multiplying x in ODE
  double m_;  // mass coefficient
  Teuchos::RCP<Teuchos::FancyOStream> out_;
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_HARMONIC_OSCILLATOR_MODEL_DECL_HPP
