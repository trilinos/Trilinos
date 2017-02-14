// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef TEMPUS_TEST_BALLPARABOLIC_MODEL_DECL_HPP
#define TEMPUS_TEST_BALLPARABOLIC_MODEL_DECL_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Tempus_Test {

/** \brief This is the "parabolic ball" model problem from Piro.
  * This is a canonical differential equation model of a ball thrown up 
  * in the air, and taking on a parabolic trajectory:
  *   \f[
  *   \ddot{x}=-1
  *   \f]
  * The initial conditions are:
  *   \f{eqnarray*}{
  *     x(0) & = & 0\\
  *     \dot{x}(0) & = & 1
  *   \f}
  * It can be shown that the exact solution to this problem is:
  *    \f[ 
  *    x(t) = t(1-0.5t)
  *    \f]
  * We consider the problem for \f$t\in [0,2]\f$ .  

  * Therefore this model has three model parameters and two initial conditions
  * which effect the exact solution as above.
  * \f[
  *   \mathbf{p}=(a,f,L)
  * \f]
  * \f[
  *   \dot{\mathbf{x}}=\mathbf{F}(\mathbf{x},t,\mathbf{p})
  * \f]
  * where
  * \f{eqnarray*}{
  *   F_{0} & = & x_{1}\\
  *   F_{1} & = & \left(\frac{f}{L}\right)^{2}(a-x_{0})
  * \f}

  * The exact sensitivities, \f$\mathbf{s}=\partial\mathbf{x}/\partial\mathbf{p}\f$,
  * for the problem are specified as
  * \f[
  *   \mathbf{s}(t)=\left[\begin{array}{cc}
  *   1 & 0\\
  *   \left(\frac{b}{L}\right)t\,\cos\left(\left(\frac{f}{L}\right)t+\phi\right) & \left(\frac{b}{L}\right)\cos\left(\left(\frac{f}{L}\right)t+\phi\right)-\frac{b\, f\, t}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)\\
  *   -\frac{b\, f\, t}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right) & -\frac{b\, f}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{2}\, t}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)
  *   \end{array}\right]
  * \f]
  * and for the default initial conditions, \f$\phi=0\f$ and \f$b=1\f$
  * \f[
  *   \mathbf{s}(t=0)=\left[\begin{array}{cc}
  *   1 & 0\\
  *   0 & \frac{b}{L}\\
  *   0 & -\frac{f}{L^{2}}
  *   \end{array}\right]
  * \f]
  * The time differentiated sensitivities, \f$\dot{\mathbf{s}}=\partial\mathbf{s}/\partial t=\partial/\partial t(\partial\mathbf{x}/\partial\mathbf{p})=\partial/\partial\mathbf{p}(\partial\mathbf{x}/\partial t)\f$
  * are
  * \f[
  *   \dot{\mathbf{s}}(t)=\left[\begin{array}{cc}
  *   0 & 0\\
  *   \left(\frac{b}{L}\right)\cos\left(\left(\frac{f}{L}\right)t+\phi\right)-\frac{b\, f\, t}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right) & -\frac{2b\, f}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)\left(\frac{b}{L}\right)-\frac{b\, f^{2}\, t}{L^{3}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)\\
  *   -\frac{b\, f}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{2}\, t}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right) & \frac{2b\, f^{2}}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{3}\, t}{L^{4}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)
  *   \end{array}\right]
  * \f]
  */

template<class Scalar>
class BallParabolicModel
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  BallParabolicModel(Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  // Exact solution
  ModelEvaluatorBase::InArgs<Scalar> getExactSolution(double t) const;

  // Exact sensitivity solution
  ModelEvaluatorBase::InArgs<Scalar> getExactSensSolution(int j, double t) const;

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

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
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs_bar,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs_bar
    ) const;
  //@}

private:
  int dim_;         ///< Number of state unknowns (2)
  int Np_;          ///< Number of parameter vectors (1)
  int np_;          ///< Number of parameters in this vector (2)
  int Ng_;          ///< Number of observation functions (1)
  int ng_;          ///< Number of elements in this observation function (1)
  bool haveIC_;     ///< false => no nominal values are provided (default=true)
  bool acceptModelParams_; ///< Changes inArgs to require parameters
  mutable bool isInitialized_;
  mutable ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  mutable ModelEvaluatorBase::OutArgs<Scalar> outArgs_;
  mutable ModelEvaluatorBase::InArgs<Scalar>  nominalValues_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > p_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > g_space_;

  // Parameters for the model:  x_0(t) = a + b*sin(f*t+phi)
  //                            x_1(t) = b*f*cos(f*t+phi)
  Scalar a_;     ///< Model parameter
  Scalar f_;     ///< Model parameter
  Scalar L_;     ///< Model parameter
  Scalar phi_;   ///< Parameter determined from the IC
  Scalar b_;     ///< Parameter determined from the IC
  Scalar t0_ic_; ///< Time value where the initial condition is specified
};


/// Non-member constructor
//Teuchos::RCP<BallParabolicModel> sineCosineModel(
//  Teuchos::RCP<Teuchos::ParameterList> pList_)
//{
//  Teuchos::RCP<BallParabolicModel> model = rcp(new BallParabolicModel(pList_));
//  return(model);
//}


} // namespace Tempus_Test
#endif // TEMPUS_TEST_BALLPARABOLIC_MODEL_DECL_HPP
