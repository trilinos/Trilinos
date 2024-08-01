//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_VANDERPOL_IMEXPart_ImplicitMODEL_DECL_HPP
#define TEMPUS_TEST_VANDERPOL_IMEXPart_ImplicitMODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief van der Pol model formulated for the partitioned IMEX-RK.
 *
 *  This is a canonical equation of a nonlinear oscillator (Hairer, Norsett,
 *  and Wanner, pp. 111-115, and Hairer and Wanner, pp. 4-5) for an electrical
 *  circuit.  In implicit ODE form, \f$ \mathcal{F}(\dot{x},x,t) = 0 \f$,
 *  the scaled problem can be written as
 *  \f{eqnarray*}{
 *    \dot{x}_0(t) - x_1(t) & = & 0 \\
 *    \dot{x}_1(t) - [(1-x_0^2)x_1-x_0]/\epsilon & = & 0
 *  \f}
 *  where the initial conditions are
 *  \f{eqnarray*}{
 *    x_0(t_0=0) & = & 2 \\
 *    x_1(t_0=0) & = & 0
 *  \f}
 *  and the initial time derivatives are
 *  \f{eqnarray*}{
 *    \dot{x}_0(t_0=0) & = & x_1(t_0=0) = 0 \\
 *    \dot{x}_1(t_0=0) & = & [(1-x_0^2)x_1-x_0]/\epsilon = -2/\epsilon
 *  \f}
 *  For a partitioned IMEX-RK time stepper, we need to rewrite this in the
 *  following form
 *  \f{eqnarray*}{
 *    M(z,t)\, \dot{z} + G(z,t) + F(z,t) & = & 0, \\
 *    \mathcal{G}(\dot{z},z,t) + F(z,t) & = & 0,
 *  \f}
 *  where \f$\mathcal{G}(\dot{z},z,t) = M(z,t)\, \dot{z} + G(z,t)\f$,
 *  \f$M(z,t)\f$ is the mass matrix, \f$F(z,t)\f$ is the operator
 *  representing the "slow" physics (and evolved explicitly), and
 *  \f$G(z,t)\f$ is the operator representing the "fast" physics.
 *  For the van der Pol problem, we can separate the terms as follows
 *  \f[
 *    z      = \left\{\begin{array}{c} y   \\ x   \end{array}\right\}
 *           = \left\{\begin{array}{c} x_0 \\ x_1 \end{array}\right\},\;\;\;
 *    F(z,t) = \left\{\begin{array}{c} F^y(x,y,t)\\F^x(x,y,t)\end{array}\right\}
 *           = \left\{\begin{array}{c} -x_1 \\ x_0/\epsilon\end{array}\right\},
 *    \mbox{ and }
 *    G(z,t) = \left\{\begin{array}{c} 0 \\ G^x(x,y,t)\end{array}\right\}
 *           = \left\{\begin{array}{c} 0 \\
 *                   -(1-x_0^2)x_1/\epsilon \end{array}\right\}
 *  \f]
 *  where \f$M(z,t)=I\f$ is the identity matrix.
 *
 *  Thus the explicit van der Pol model (VanDerPol_IMEX_ExplicitModel)
 *  formulated for the partitioned IMEX-RK is
 *  \f{eqnarray*}{
 *    F^y(x,y,t) & = & \dot{x}_0(t) - x_1(t) = 0 \\
 *    F^x(x,y,t) & = & \dot{x}_1(t) + x_0/\epsilon = 0
 *  \f}
 *  and the implicit van der Pol model (VanDerPol_IMEXPart_ImplicitModel)
 *  formulated for the partitioned IMEX-RK is
 *  \f[
 *    G^x(x,y,t) = \dot{x}_1(t) - (1-x_0^2)x_1/\epsilon = 0
 *  \f]
 *  Noting that \f$G^y(x,y,t) = \dot{x}_0(t) = 0\f$ is not needed.
 *
 *  Recalling the defintion of the iteration matrix, \f$W\f$,
 *  \f[
 *    W_{ij} \equiv \frac{d\mathcal{G}^x_i}{dx_j} =
 *      \alpha \frac{\partial\mathcal{G}^x_i}{\partial \dot{x}_j}
 *    + \beta \frac{\partial\mathcal{G}^x_i}{\partial x_j}
 *  \f]
 *  where
 *  \f[
 *    \alpha = \left\{
 *      \begin{array}{cl}
 *        \frac{\partial\dot{x}_i}{\partial x_j} & \mbox{ if } i = j \\
 *        0 & \mbox{ if } i \neq j
 *      \end{array} \right.
 *    \;\;\;\; \mbox{ and } \;\;\;\;
 *    \beta = 1
 *  \f]
 *  we can write for the implicit van der Pol model
 *  (VanDerPol_IMEXPart_ImplicitModel)
 *  \f{eqnarray*}{
 *    W_{00} = \alpha \frac{\partial\mathcal{G}^x_0}{\partial \dot{x}_1}
 *            + \beta \frac{\partial\mathcal{G}^x_0}{\partial x_1}
 *         & = & \alpha + \beta (x^2_0 - 1)/\epsilon \\
 *  \f}
 */

template <class Scalar>
class VanDerPol_IMEXPart_ImplicitModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase {
 public:
  // Constructor
  VanDerPol_IMEXPart_ImplicitModel(
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
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar) const;
  //@}

  int dim_;                ///< Number of state unknowns (1)
  int Np_;                 ///< Number of parameter vectors (1)
  int np_;                 ///< Number of parameters in this vector (1)
  int Ng_;                 ///< Number of observation functions (0)
  int ng_;                 ///< Number of elements in this observation function (0)
  bool haveIC_;            ///< false => no nominal values are provided (default=true)
  bool useDfDpAsTangent_;  ///< Treat DfDp OutArg as tangent (df/dx*dx/dp+df/dp)
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > y_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > dxdp_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > dydp_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  // Parameters for the model:
  Scalar epsilon_;  ///< This is a model parameter
  Scalar t0_ic_;    ///< initial time
  Scalar x0_ic_;    ///< initial condition for x0
  Scalar x1_ic_;    ///< initial condition for x1
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_VANDERPOL_IMEXPart_ImplicitMODEL_DECL_HPP
