// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TEMPUS_SINCOS_MODEL_DECL_HPP
#define ROL_TEMPUS_SINCOS_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

/** \brief Sine-Cosine model problem from Rythmos.
  * This is a canonical Sine-Cosine differential equation
  *   \f[
  *   \mathbf{\ddot{x}}=-\mathbf{x}
  *   \f]
  * with a few enhancements. We start with the exact solution to the
  * differential equation
  *   \f{eqnarray*}{
  *     x_{0}(t) & = & a+b*\sin(w*t+\phi)\\
  *     x_{1}(t) & = & b*w*\cos(w*t+\phi)
  *   \f}
  * then the form of the model is
  * \f{eqnarray*}{
  *   \frac{d}{dt}x_{0}(t) & = & x_{1}(t)\\
  *   \frac{d}{dt}x_{1}(t) & = & w^{2}(a-x_{0}(t))
  * \f}
  * where the default parameter values are \f$a=0\f$, and \f$w=1\f$,
  * and the initial conditions
  * \f{eqnarray*}{
  *   x_{0}(t_{0}=0) & = & \gamma_{0}[=0]\\
  *   x_{1}(t_{0}=0) & = & \gamma_{1}[=1]
  * \f}
  * determine the remaining coefficients
  * \f{eqnarray*}{
  *   \phi & = & \arctan((w/\gamma_{1})*(\gamma_{0}-a))-w*t_{0}[=0]\\
  *   b & = & \gamma_{1}/(w*cos(w*t_{0}+\phi))[=1]
  * \f}

  * Therefore this model has two model parameters and two initial conditions
  * which effect the exact solution as above.
  * \f[
  *   \mathbf{p}=(a,w)
  * \f]
  * \f[
  *   \dot{\mathbf{x}}=\mathbf{F}(\mathbf{x},t,\mathbf{p})
  * \f]
  * where
  * \f{eqnarray*}{
  *   F_{0} & = & x_{1}\\
  *   F_{1} & = & w^{2}(a-x_{0})
  * \f}
  */

template<class Scalar>
class SinCosModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  SinCosModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  // Exact sensitivity solution
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getExactSensSolution(int j, double t) const;

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;

  //@}

  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

private:

  void setupInOutArgs_() const;

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
    ) const;
  //@}

  void calculateCoeffFromIC_();

private:
  int dim_;         ///< Number of state unknowns (2)
  int Np_;          ///< Number of parameter vectors (1)
  int np_;          ///< Number of parameters in this vector (2)
  int Ng_;          ///< Number of observation functions (1)
  int ng_;          ///< Number of elements in this observation function (1)
  bool haveIC_;     ///< false => no nominal values are provided (default=true)
  bool acceptModelParams_; ///< Changes inArgs to require parameters
  bool useDfDpAsTangent_; ///< Treat DfDp OutArg as tangent (df/dx*dx/dp+df/dp)
  mutable bool isInitialized_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable Thyra::ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > DxDp_space_;

  // Parameters for the model:  x_0(t) = a + b*sin(w*t+phi)
  //                            x_1(t) = b*w*cos(w*t+phi)
  Scalar a_;     ///< Model parameter
  Scalar w_;     ///< Model parameter
  Scalar phi_;   ///< Parameter determined from the IC
  Scalar b_;     ///< Parameter determined from the IC
  Scalar t0_ic_; ///< Time value where the initial condition is specified
  Scalar x0_ic_; ///< Initial condition for x0
  Scalar x1_ic_; ///< Initial condition for x1
};

#endif //  ROL_TEMPUS_SINCOS_MODEL_DECL_HPP

#include "example_sincos_impl.hpp"
