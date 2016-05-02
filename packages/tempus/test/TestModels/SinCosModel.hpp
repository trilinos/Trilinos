#ifndef TEMPUS_TEST_SINCOS_MODEL_HPP
#define TEMPUS_TEST_SINCOS_MODEL_HPP

#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_ParameterList.hpp"

using Thyra::ModelEvaluatorBase;

namespace Tempus_Test {

/** \brief Sine-Cosine model problem from Rythmos.
  * This is a canonical Sine-Cosine differential equation
  *   \f[
  *   \mathbf{\ddot{x}}=-\mathbf{x}
  *   \f]
  * with a few enhancements. We start with the exact solution to the
  * differential equation
  *   \f{eqnarray*}{
  *     x_{0}(t) & = & a+b*\sin((f/L)*t+\phi)\\
  *     x_{1}(t) & = & b*(f/L)*\cos((f/L)*t+\phi)
  *   \f}
  * then the form of the model is
  * \f{eqnarray*}{
  *   \frac{d}{dt}x_{0}(t) & = & x_{1}(t)\\
  *   \frac{d}{dt}x_{1}(t) & = & \left(\frac{f}{L}\right)^{2}(a-x_{0}(t))
  * \f}
  * where the default parameter values are $a=0$, $f=1$, and $L=1$,
  * and the initial conditions
  * \f{eqnarray*}{
  *   x_{0}(t_{0}=0) & = & \gamma_{0}[=0]\\
  *   x_{1}(t_{0}=0) & = & \gamma_{1}[=1]
  * \f}
  * determine the remaining coefficients
  * \f{eqnarray*}{
  *   \phi & = & \arctan(((f/L)/\gamma_{1})*(\gamma_{0}-a))-(f/L)*t_{0}[=0]\\
  *   b & = & \gamma_{1}/((f/L)*cos((f/L)*t_{0}+\phi))[=1]
  * \f}

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

  * The exact sensitivities, $\mathbf{s}=\partial\mathbf{x}/\partial\mathbf{p}$,
  * for the problem are specified as
  * \f[
  *   \mathbf{s}(t)=\left[\begin{array}{cc}
  *   1 & 0\\
  *   \left(\frac{b}{L}\right)t\,\cos\left(\left(\frac{f}{L}\right)t+\phi\right) & \left(\frac{b}{L}\right)\cos\left(\left(\frac{f}{L}\right)t+\phi\right)-\frac{b\, f\, t}{L^{2}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)\\
  *   -\frac{b\, f\, t}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right) & -\frac{b\, f}{L^{2}}\cos\left(\left(\frac{f}{L}\right)t+\phi\right)+\frac{b\, f^{2}\, t}{L^{3}}\sin\left(\left(\frac{f}{L}\right)t+\phi\right)
  *   \end{array}\right]
  * \f]
  * and for the default initial conditions, $\phi=0$ and $b=1$
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

class SinCosModel
  : public Thyra::StateFuncModelEvaluatorBase<double>,
    public Teuchos::ParameterListAcceptorDefaultBase
{
  public:

  // Constructor
  SinCosModel(Teuchos::RCP<Teuchos::ParameterList> pList);

  // Exact solution
  ModelEvaluatorBase::InArgs<double> getExactSolution(double t) const;

  // Exact sensitivity solution
  ModelEvaluatorBase::InArgs<double> getExactSensSolution(int j, double t) const;

  // Set explicit/implicit flag
  void setImplicitFlag(bool implicit);

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const;
  RCP<const Thyra::VectorSpaceBase<double> > get_f_space() const;
  ModelEvaluatorBase::InArgs<double> getNominalValues() const;
  RCP<Thyra::LinearOpWithSolveBase<double> > create_W() const;
  RCP<Thyra::LinearOpBase<double> > create_W_op() const;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > get_W_factory() const;
  ModelEvaluatorBase::InArgs<double> createInArgs() const;

  RCP<const Thyra::VectorSpaceBase<double> > get_p_space(int l) const;
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  RCP<const Thyra::VectorSpaceBase<double> > get_g_space(int j) const;

  //@}

  /** \name Public functions overridden from ParameterListAcceptor. */
  //@{
  void setParameterList(RCP<ParameterList> const& paramList);
  RCP<const ParameterList> getValidParameters() const;
  //@}

private:

  void setupInOutArgs() const;

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{
  ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const;
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<double> &inArgs_bar,
    const ModelEvaluatorBase::OutArgs<double> &outArgs_bar
    ) const;
  //@}

  void calculateCoeffFromIC();

private:
  int dim;         ///< Number of state unknowns (2)
  int Np;          ///< Number of parameter vectors (1)
  int np;          ///< Number of parameters in this vector (2)
  int Ng;          ///< Number of observation functions (1)
  int ng;          ///< Number of elements in this observation function (1)
  bool isImplicit; ///< false => xdot = f(x,t)   W = beta*df/dx; true =>  F(xdot,x,t) = 0 W = alpha*dF/dxdot + beta*dF/dx
  bool haveIC;     ///< false => no nominal values are provided (default=true)
  bool acceptModelParams; ///< Changes inArgs to require parameters
  mutable bool isInitialized;
  mutable ModelEvaluatorBase::InArgs<double>  inArgs;
  mutable ModelEvaluatorBase::OutArgs<double> outArgs;
  mutable ModelEvaluatorBase::InArgs<double>  nominalValues;
  RCP<const Thyra::VectorSpaceBase<double> >  x_space;
  RCP<const Thyra::VectorSpaceBase<double> >  f_space;
  RCP<const Thyra::VectorSpaceBase<double> >  p_space;
  RCP<const Thyra::VectorSpaceBase<double> >  g_space;

  // Parameters for the model:  x_0(t) = a + b*sin(f*t+phi)
  //                            x_1(t) = b*f*cos(f*t+phi)
  double a_;    ///< Model parameter
  double f_;    ///< Model parameter
  double L_;    ///< Model parameter
  double phi_;  ///< Parameter determined from the IC
  double b_;    ///< Parameter determined from the IC
  double t0_ic; ///< Time value where the initial condition is specified
  double x0_ic; ///< Initial condition for x0
  double x1_ic; ///< Initial condition for x1
};


// Non-member constructor
RCP<SinCosModel> sinCosModel(Teuchos::RCP<Teuchos::ParameterList> pList);


} // namespace Tempus_Test
#endif // TEMPUS_TEST_SINCOS_MODEL_HPP
