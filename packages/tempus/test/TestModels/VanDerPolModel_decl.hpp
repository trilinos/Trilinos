//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_VANDERPOL_MODEL_DECL_HPP
#define TEMPUS_TEST_VANDERPOL_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp"               // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp"  // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Tempus_Test {

/** \brief van der Pol model problem for nonlinear electrical circuit.
 *
 * This is a canonical equation of a nonlinear oscillator (Hairer, Norsett,
 * and Wanner, pp. 111-115, and Hairer and Wanner, pp. 4-5) for an electrical
 * circuit.  In implicit ODE form, \f$ \mathcal{F}(\dot{x},x,t) = 0 \f$,
 * the scaled problem can be written as
 * \f{eqnarray*}{
 *   \mathcal{F}_0 & = & \dot{x}_0(t) - x_1(t) = 0 \\
 *   \mathcal{F}_1 & = & \dot{x}_1(t) - [(1-x_0^2)x_1-x_0]/\epsilon = 0
 * \f}
 * where the initial conditions are
 * \f{eqnarray*}{
 *   x_0(t_0=0) & = & 2 \\
 *   x_1(t_0=0) & = & 0
 * \f}
 * and the initial time derivatives are
 * \f{eqnarray*}{
 *   \dot{x}_0(t_0=0) & = & x_1(t_0=0) = 0 \\
 *   \dot{x}_1(t_0=0) & = & [(1-x_0^2)x_1-x_0]/\epsilon = -2/\epsilon
 * \f}
 * Hairer and Wanner suggest the output times of \f$t = 1,2,3,4,...,11\f$,
 * and \f$\epsilon = 10^{-6}\f$ to make the problem very stiff.
 * For \f$\epsilon = 0\f$, the solution becomes
 * \f{eqnarray*}{
 *   \ln \left|x_0\right| - \frac{x_0^2}{2} & = & t + C \\
 *   x_1 & = & \frac{x_0}{1-x_0^2}
 * \f}
 * where \f$C =\ln \left|x_0(t=0)\right| - \frac{x_0^2(t=0)}{2} =-1.306853.\f$
 *
 * The components of iteration matrix, \f$W\f$, are defined to be
 * \f[
 *   W_{ij} \equiv \frac{d\mathcal{F}_i}{dx_j} = \frac{d}{dx_j}
 *          \mathcal{F}_i (\dot{x}_i, x_0, \ldots, x_k, \ldots, x_K, t)
 * \f]
 * (not using Einstein summation).  Using the chain rule, we can write
 * \f[
 *   \frac{d\mathcal{F}_i}{dx_j} =
 *    \frac{\partial\dot{x}_i}{\partial x_j}
 *    \frac{\partial\mathcal{F}_i}{\partial \dot{x}_i}
 *    + \sum_{k=0}^K \frac{\partial x_k}{\partial x_j}
 *      \frac{\partial\mathcal{F}_i}{\partial x_k}
 *    + \frac{\partial t}{\partial x_j}
 *    \frac{\partial\mathcal{F}_i}{\partial t}
 * \f]
 * but noting that \f$\partial t/\partial x_j = 0\f$ and
 * \f[
 *    \frac{\partial x_k}{\partial x_j} = \left\{
 *      \begin{array}{c}
 *        1 \mbox{ if } j = k \\
 *        0 \mbox{ if } j \neq k
 *      \end{array}
 *    \right.
 * \f]
 * we can write
 * \f[
 *   \frac{d\mathcal{F}_i}{dx_j} =
 *     \alpha \frac{\partial\mathcal{F}_i}{\partial \dot{x}_j}
 *   + \beta \frac{\partial\mathcal{F}_i}{\partial x_j}
 * \f]
 * where
 * \f[
 *   \alpha = \left\{
 *     \begin{array}{cl}
 *       \frac{\partial\dot{x}_i}{\partial x_j} & \mbox{ if } i = j \\
 *       0 & \mbox{ if } i \neq j
 *     \end{array} \right.
 *   \;\;\;\; \mbox{ and } \;\;\;\;
 *   \beta = \left\{
 *   \begin{array}{cl}
 *     \frac{\partial x_k}{\partial x_j} = 1 & \mbox{ if } j = k \\
 *     0 & \mbox{ if } j \neq k
 *   \end{array} \right.
 * \f]
 * Thus for the van der Pol problem, we have
 * \f{eqnarray*}{
 *   W_{00} = \alpha \frac{\partial\mathcal{F}_0}{\partial \dot{x}_0}
 *           + \beta \frac{\partial\mathcal{F}_0}{\partial x_0}
 *        & = & \alpha \\
 *   W_{01} = \alpha \frac{\partial\mathcal{F}_0}{\partial \dot{x}_1}
 *           + \beta \frac{\partial\mathcal{F}_0}{\partial x_1}
 *        & = & -\beta \\
 *   W_{10} = \alpha \frac{\partial\mathcal{F}_1}{\partial \dot{x}_0}
 *           + \beta \frac{\partial\mathcal{F}_1}{\partial x_0}
 *        & = & \beta (2 x_0 x_1 + 1)/\epsilon \\
 *   W_{11} = \alpha \frac{\partial\mathcal{F}_1}{\partial \dot{x}_1}
 *           + \beta \frac{\partial\mathcal{F}_1}{\partial x_1}
 *        & = & \alpha + \beta (x^2_0 - 1)/\epsilon \\
 * \f}
 */

template <class Scalar>
class VanDerPolModel : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
                       public Teuchos::ParameterListAcceptorDefaultBase {
 public:
  // Constructor
  VanDerPolModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  // Exact sensitivity solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSensSolution(
      int j, double t) const;

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

  int dim_;                 ///< Number of state unknowns (2)
  int Np_;                  ///< Number of parameter vectors (1)
  int np_;                  ///< Number of parameters in this vector (1)
  int Ng_;                  ///< Number of observation functions (0)
  int ng_;                  ///< Number of elements in this observation function (0)
  bool haveIC_;             ///< false => no nominal values are provided (default=true)
  bool acceptModelParams_;  ///< Changes inArgs to require parameters
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  // Parameters for the model:
  Scalar epsilon_;  ///< This is a model parameter
  Scalar t0_ic_;    ///< initial time
  Scalar x0_ic_;    ///< initial condition for x0
  Scalar x1_ic_;    ///< initial condition for x1
};

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_VANDERPOL_MODEL_DECL_HPP
