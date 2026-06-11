//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_LOTKAVOLTERRA_MODEL_DECL_HPP
#define TEMPUS_TEST_LOTKAVOLTERRA_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief Lotka-Volterra predator-prey model with optional time-dependent
 *  forcing on the prey equation.
 *
 * This is the Lotka-Volterra system describing the interaction between a
 * prey population \f$x_0\f$ and a predator population \f$x_1\f$, augmented
 * with a sinusoidal forcing term on the prey equation:
 * \f{eqnarray*}{
 *   \frac{d}{dt}x_0(t) & = & \alpha\, x_0 - \beta\, x_0\, x_1
 *                            - k\,\sin(t) \\
 *   \frac{d}{dt}x_1(t) & = & \delta\, x_0\, x_1 - \gamma\, x_1
 * \f}
 * The default parameter values are
 * \f$\alpha=1.5\f$, \f$\beta=1.0\f$, \f$\delta=1.0\f$, \f$\gamma=3.0\f$,
 * \f$k=0\f$ (unforced), and the default initial conditions are
 * \f{eqnarray*}{
 *   x_0(t_0=0) & = & 10.0 \\
 *   x_1(t_0=0) & = & 5.0
 * \f}
 *
 * When \f$k=0\f$ the system reduces to the classical autonomous
 * Lotka-Volterra equations, which possess a conserved quantity (first
 * integral):
 * \f[
 *   H(x_0, x_1) = \delta\, x_0 - \gamma\, \ln(x_0)
 *                + \beta\, x_1  - \alpha\, \ln(x_1).
 * \f]
 * For \f$k \neq 0\f$ this quantity is no longer conserved.
 *
 * The equilibrium point of the unforced system is
 * \f$\left(x_0^*, x_1^*\right) = \left(\gamma/\delta,\; \alpha/\beta\right)\f$.
 *
 * In implicit ODE form \f$\mathcal{F}(\dot{x}, x, t) = 0\f$:
 * \f{eqnarray*}{
 *   \mathcal{F}_0 & = & \dot{x}_0 - \alpha\,x_0 + \beta\,x_0\,x_1
 *                       + k\,\sin(t) = 0 \\
 *   \mathcal{F}_1 & = & \dot{x}_1 - \delta\,x_0\,x_1 + \gamma\,x_1 = 0
 * \f}
 *
 * Because the forcing is state-independent, the iteration matrix
 * \f$W = \alpha_c\,\partial\mathcal{F}/\partial\dot{x}
 * + \beta_c\,\partial\mathcal{F}/\partial x\f$ is unchanged from the
 * unforced case:
 * \f{eqnarray*}{
 *   W_{00} & = & \alpha_c + \beta_c\,(\beta\,x_1 - \alpha) \\
 *   W_{01} & = & \beta_c\,\beta\,x_0 \\
 *   W_{10} & = & \beta_c\,(-\delta\,x_1) \\
 *   W_{11} & = & \alpha_c + \beta_c\,(\delta\,x_0 - \gamma)
 * \f}
 * (The subscript \f$c\f$ on \f$\alpha_c, \beta_c\f$ distinguishes the
 * Tempus time-integration scalars from the Lotka-Volterra model coefficients
 * \f$\alpha, \beta\f$.)
 *
 * For the explicit ODE (\f$\dot{x}\f$ is null), the Jacobian reduces to
 * \f$\beta_c \cdot \partial f/\partial x\f$:
 * \f{eqnarray*}{
 *   W_{00} & = & \beta_c\,(\alpha - \beta\,x_1) \\
 *   W_{01} & = & \beta_c\,(-\beta\,x_0) \\
 *   W_{10} & = & \beta_c\,\delta\,x_1 \\
 *   W_{11} & = & \beta_c\,(\delta\,x_0 - \gamma)
 * \f}
 */
template <class Scalar>
class LotkaVolterraModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase {
 public:
  // Constructor
  LotkaVolterraModel(
      Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

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
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;
  //@}

 protected:
  int  dim_;                ///< Number of state unknowns (2)
  int  Np_;                 ///< Number of parameter vectors (0)
  int  Ng_;                 ///< Number of observation functions (0)
  bool haveIC_;             ///< false => no nominal values provided (default=true)
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;

  // Lotka-Volterra model parameters
  Scalar alpha_;   ///< Prey   growth rate
  Scalar beta_;    ///< Predation rate
  Scalar delta_;   ///< Predator growth rate per prey killed
  Scalar gamma_;   ///< Predator death rate
  Scalar k_;       ///< Forcing amplitude on prey
  Scalar x0_ic_;   ///< Initial condition for prey  x0
  Scalar y0_ic_;   ///< Initial condition for predator x1
  Scalar t0_ic_;   ///< Initial time t0
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_LOTKAVOLTERRA_MODEL_DECL_HPP
