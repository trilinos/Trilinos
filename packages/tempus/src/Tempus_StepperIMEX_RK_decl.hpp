// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperIMEX_RK_decl_hpp
#define Tempus_StepperIMEX_RK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"
#include "Tempus_StepperIMEX_RKObserver.hpp"


namespace Tempus {

/** \brief Implicit-Explicit Runge-Kutta (IMEX-RK) time stepper.
 *
 *  For the implicit ODE system, \f$ \mathcal{F}(\dot{x},x,t) = 0 \f$,
 *  we need to specialize this in order to separate the explicit, implicit,
 *  and temporal terms for the IMEX-RK time stepper,
 *  \f{eqnarray*}{
 *    M(x,t)\, \dot{x}(x,t) + G(x,t) + F(x,t) & = & 0, \\
 *    \mathcal{G}(\dot{x},x,t) + F(x,t) & = & 0,
 *  \f}
 *  where \f$\mathcal{G}(\dot{x},x,t) = M(x,t)\, \dot{x} + G(x,t)\f$,
 *  \f$M(x,t)\f$ is the mass matrix, \f$F(x,t)\f$ is the operator
 *  representing the "slow" physics (and is evolved explicitly), and
 *  \f$G(x,t)\f$ is the operator representing the "fast" physics (and is
 *  evolved implicitly).  Additionally, we assume that the mass matrix is
 *  invertible, so that
 *  \f[
 *    \dot{x}(x,t) + g(x,t) + f(x,t) = 0
 *  \f]
 *  where \f$f(x,t) = M(x,t)^{-1}\, F(x,t)\f$, and
 *  \f$g(x,t) = M(x,t)^{-1}\, G(x,t)\f$. Using Butcher tableaus for the
 *  explicit terms
 *  \f[ \begin{array}{c|c}
 *    \hat{c} & \hat{A} \\ \hline
 *            & \hat{b}^T
 *  \end{array}
 *  \;\;\;\; \mbox{ and for implicit terms } \;\;\;\;
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}, \f]
 *  the basic IMEX-RK method for \f$s\f$-stages can be written as
 *  \f[ \begin{array}{rcll}
 *    X_i & = & x_{n-1}
 *     - \Delta t \sum_{j=1}^{i-1} \hat{a}_{ij}\, f(X_j,\hat{t}_j)
 *     - \Delta t \sum_{j=1}^i           a_{ij}\, g(X_j,t_j)
 *            & \mbox{for } i=1\ldots s, \\
 *    x_n & = & x_{n-1}
 *     - \Delta t \sum_{i=1}^s \hat{b}_{i}\, f(X_i,\hat{t}_i)
 *     - \Delta t \sum_{i=1}^s       b_{i}\, g(X_i,t_i) &
 *  \end{array} \f]
 *  where \f$\hat{t}_i = t_{n-1}+\hat{c}_i\Delta t\f$ and
 *  \f$t_i = t_{n-1}+c_i\Delta t\f$.  For iterative solvers, it is useful to
 *  write the stage solutions as
 *  \f[
 *    X_i = \tilde{X} - a_{ii} \Delta t\, g(X_i,t_i)
 *  \f]
 *  where
 *  \f[
 *    \tilde{X} = x_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *        \left(\hat{a}_{ij}\, f(X_j,\hat{t}_j) + a_{ij}\, g(X_j,t_j)\right)
 *  \f]
 *  Rearranging to solve for the implicit term
 *  \f[
 *    g(X_i,t_i) = - \frac{X_i - \tilde{X}}{a_{ii} \Delta t}
 *  \f]
 *  We can use this to determine the time derivative at each stage for the
 *  implicit solve.
 *  \f[
 *    \dot{X}_i(X_i,t_i) + g(X_i,t_i) + f(X_i,t_i) = 0
 *  \f]
 *  Note that the explicit term, \f$f(X_i,t_i)\f$, is evaluated at the implicit
 *  stage time, \f$t_i\f$.
 *  We can form the time derivative
 *  \f{eqnarray*}{
 *    \dot{X}(X_i,t_i) & = & - g(X_i,t_i) - f(X_i,t_i) \\
 *    \dot{X}(X_i,t_i) & = &
 *      \frac{X_i - \tilde{X}}{a_{ii} \Delta t} - f(X_i,t_i) \\
 *    \dot{X}(X_i,t_i) & = &
 *      \frac{X_i - \tilde{X}}{a_{ii} \Delta t} -M(X_i, t_i)^{-1}\, F(X_i,t_i)\\
 *  \f}
 *  Returning to the governing equation
 *  \f{eqnarray*}{
 *    M(X_i,t_i)\, \dot{X}(X_i,t_i) + G(X_i,t_i) + F(X_i,t_i) & = & 0 \\
 *    M(X_i,t_i)\, \left[
 *      \frac{X_i - \tilde{X}}{a_{ii} \Delta t} - M(X_i, t_i)^{-1}\, F(X_i,t_i)
 *    \right] + G(X_i,t_i) + F(X_i,t_i) & = & 0 \\
 *      M(X_i,t_i)\, \left[ \frac{X_i - \tilde{X}}{a_{ii} \Delta t} \right]
 *    + G(X_i,t_i) & = & 0 \\
 *  \f}
 *  Recall \f$\mathcal{G}(\dot{x},x,t) = M(x,t)\, \dot{x} + G(x,t)\f$ and if we
 *  define a pseudo time derivative,
 *  \f[
 *    \tilde{\dot{X}} = \frac{X_i - \tilde{X}}{a_{ii} \Delta t} = - g(X_i,t_i),
 *  \f]
 *  we can write
 *  \f[
 *    \mathcal{G}(\tilde{\dot{X}},X_i,t_i) =
 *       M(X_i,t_i)\, \tilde{\dot{X}} + G(X_i,t_i) = 0
 *  \f]
 *
 *  For the case when \f$a_{ii}=0\f$, we need the time derivative for the
 *  plain explicit case.  Thus the stage solution is
 *  \f[
 *     X_i = x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *            \hat{a}_{ij}\, f_j + a_{ij}\, g_j \right) = \tilde{X}
 *  \f]
 *  and we can simply evaluate
 *  \f{eqnarray*}{
 *     f_i & = & M(X_i,\hat{t}_i)^{-1}\, F(X_i,\hat{t}_i) \\
 *     g_i & = & M(X_i,      t_i)^{-1}\, G(X_i,      t_i)
 *  \f}
 *  We can then form the time derivative as
 *  \f[
 *    \dot{X}_i(X_i,t_i) = - g(X_i,t_i) - f(X_i,t_i)
 *  \f]
 *  but again note that the explicit term, \f$f(X_i,t_i)\f$, is evaluated
 *  at the implicit stage time, \f$t_i\f$.
 *
 *  <b> IMEX-RK Algorithm </b>
 *
 *  The single-timestep algorithm for IMEX-RK using the real time derivative,
 *  \f$\dot{X}(X_i,t_i)\f$, is
 *   - \f$X_1 \leftarrow x_{n-1}\f$
 *   - for \f$i = 1 \ldots s\f$ do
 *     - \f$\tilde{X} \leftarrow x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *            \hat{a}_{ij}\, f_j + a_{ij}\, g_j \right) \f$
 *     - if \f$a_{ii} = 0\f$
 *       - \f$X_i \leftarrow \tilde{X}\f$
 *       - \f$g_i \leftarrow M(X_i,      t_i)^{-1}\, G(X_i,      t_i)\f$
 *     - else
 *       - Define \f$\dot{X}(X_i,t_i)
 *            = \frac{X_i-\tilde{X}}{a_{ii} \Delta t}
 *                     - M(X_i,t_i)^{-1}\, F(X_i,t_i) \f$
 *       - Solve \f$\mathcal{G}\left(\dot{X}(X_i,t_i),X_i,t_i\right)
 *                  + F(X_i,t_i) = 0\f$ for \f$X_i\f$
 *       - \f$g_i \leftarrow - \frac{X_i-\tilde{X}}{a_{ii} \Delta t}\f$
 *     - \f$f_i \leftarrow M(X_i,\hat{t}_i)^{-1}\, F(X_i,\hat{t}_i)\f$
 *     - \f$\dot{X}_i = - g_i - M(X_i,t_i)^{-1}\,F(X_i,t_i)\f$ [Optional]
 *   - end for
 *   - \f$x_n \leftarrow x_{n-1} - \Delta t\,\sum_{i=1}^{s}\hat{b}_i\,f_i
 *                               - \Delta t\,\sum_{i=1}^{s}      b_i\,g_i\f$
 *   - Solve \f$M(x_n)\, \dot{x}_n + F(x_n,t_n) + G(x_n,t_n) = 0\f$
 *       for \f$\dot{x}_n\f$ [Optional]
 *
 *  The single-timestep algorithm for IMEX-RK using the pseudo time derivative,
 *  \f$\tilde{\dot{X}}\f$, is (which is currently implemented)
 *   - \f$X_1 \leftarrow x_{n-1}\f$
 *   - for \f$i = 1 \ldots s\f$ do
 *     - \f$\tilde{X} \leftarrow x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *            \hat{a}_{ij}\, f_j + a_{ij}\, g_j \right) \f$
 *     - if \f$a_{ii} = 0\f$
 *       - \f$X_i \leftarrow \tilde{X}\f$
 *       - \f$g_i \leftarrow M(X_i,      t_i)^{-1}\, G(X_i,      t_i)\f$
 *     - else
 *       - Define \f$\tilde{\dot{X}}
 *            = \frac{X_i-\tilde{X}}{a_{ii} \Delta t} \f$
 *       - Solve \f$\mathcal{G}\left(\tilde{\dot{X}},X_i,t_i\right) = 0\f$
 *           for \f$X_i\f$
 *       - \f$g_i \leftarrow - \tilde{\dot{X}}\f$
 *     - \f$f_i \leftarrow M(X_i,\hat{t}_i)^{-1}\, F(X_i,\hat{t}_i)\f$
 *     - \f$\dot{X}_i = - g_i - M(X_i,t_i)^{-1}\,F(X_i,t_i)\f$ [Optional]
 *   - end for
 *   - \f$x_n \leftarrow x_{n-1} - \Delta t\,\sum_{i=1}^{s}\hat{b}_i\,f_i
 *                               - \Delta t\,\sum_{i=1}^{s}      b_i\,g_i\f$
 *   - Solve \f$M(x_n)\, \dot{x}_n + F(x_n,t_n) + G(x_n,t_n) = 0\f$
 *       for \f$\dot{x}_n\f$ [Optional]
 *
 *  The following table contains the pre-coded IMEX-RK tableaus.
 *  <table>
 *  <caption id="multi_row">IMEX-RK Tableaus</caption>
 *  <tr><th> Name  <th> Order <th> Implicit Tableau <th> Explicit Tableau
 *  <tr><td> IMEX RK 1st order  <td> 1st
 *      <td> \f[ \begin{array}{c|cc}
 *             0 & 0 & 0 \\
 *             1 & 0 & 1 \\ \hline
 *               & 0 & 1
 *           \end{array} \f]
 *      <td> \f[ \begin{array}{c|cc}
 *             0 & 0 & 0 \\
 *             1 & 1 & 0 \\ \hline
 *               & 1 & 0
 *           \end{array} \f]
 *  <tr><td> IMEX RK SSP2\n \f$\gamma = 1-1/\sqrt{2}\f$ <td> 2nd
 *      <td> \f[ \begin{array}{c|cc}
 *             \gamma   & \gamma & 0 \\
 *             1-\gamma & 1-2\gamma & \gamma \\ \hline
 *                      & 1/2       & 1/2
 *           \end{array} \f]
 *      <td> \f[ \begin{array}{c|cc}
 *             0 & 0   & 0 \\
 *             1 & 1   & 0 \\ \hline
 *               & 1/2 & 1/2
 *           \end{array} \f]
 *  <tr><td> IMEX RK ARS 233\n \f$\gamma = (3+\sqrt{3})/6\f$ <td> 3rd
 *      <td> \f[ \begin{array}{c|ccc}
 *             0        & 0      & 0         & 0      \\
 *             \gamma   & 0      & \gamma    & 0      \\
 *             1-\gamma & 0      & 1-2\gamma & \gamma \\ \hline
 *                      & 0      & 1/2       & 1/2
 *           \end{array} \f]
 *      <td> \f[ \begin{array}{c|ccc}
 *             0        & 0        & 0         & 0 \\
 *             \gamma   & \gamma   & 0         & 0 \\
 *             1-\gamma & \gamma-1 & 2-2\gamma & 0 \\ \hline
 *                      & 0        & 1/2       & 1/2
 *           \end{array} \f]
 *  </table>
 *
 *  #### References
 *  -# Ascher, Ruuth, Spiteri, "Implicit-explicit Runge-Kutta methods for
 *     time-dependent partial differential equations", Applied Numerical
 *     Mathematics 25 (1997) 151-167.
 *  -# Cyr, "IMEX Lagrangian Methods", SAND2015-3745C.
 *  -# Shadid, Cyr, Pawlowski, Widley, Scovazzi, Zeng, Phillips, Conde,
 *     Chuadhry, Hensinger, Fischer, Robinson, Rider, Niederhaus, Sanchez,
 *     "Towards an IMEX Monolithic ALE Method with Integrated UQ for
 *     Multiphysics Shock-hydro", SAND2016-11353, 2016, pp. 21-28.
 */
template<class Scalar>
class StepperIMEX_RK : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor to use default Stepper parameters.
  StepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    std::string stepperType = "IMEX RK SSP2");

  /// Constructor to specialize Stepper parameters.
  StepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList);

  /// Constructor for StepperFactory.
  StepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& models,
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> pList);

  /// \name Basic stepper methods
  //@{
    /// Set both the explicit and implicit tableau from ParameterList
    virtual void setTableaus(Teuchos::RCP<Teuchos::ParameterList> pList,
                             std::string stepperType = "");

    /// Set the explicit tableau from ParameterList
    virtual void setExplicitTableau(std::string stepperType,
                                    Teuchos::RCP<Teuchos::ParameterList> pList);

    /// Set the explicit tableau from tableau
    virtual void setExplicitTableau(
      Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau);

    /// Set the implicit tableau from ParameterList
    virtual void setImplicitTableau(std::string stepperType,
                                    Teuchos::RCP<Teuchos::ParameterList> pList);

    /// Set the implicit tableau from tableau
    virtual void setImplicitTableau(
      Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau);

    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel()
     { return this->wrapperModel_; }

    virtual void setModelPair(
      const Teuchos::RCP<WrapperModelEvaluatorPairIMEX_Basic<Scalar> >& mePair);

    virtual void setModelPair(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel);

    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Pass initial guess to Newton solver (only relevant for implicit solvers)
    virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
       {initial_guess_ = initial_guess;}

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> >getDefaultStepperState();
    virtual Scalar getOrder()const { return order_; }
    virtual Scalar getOrderMin()const { return order_; }
    virtual Scalar getOrderMax()const { return order_; }

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
  //@}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> getDefaultParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual std::string description() const;
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  void evalImplicitModelExplicitly(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
    Scalar time, Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & G) const;

  void evalExplicitModel(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
    Scalar time, Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const;

private:

  /// Default Constructor -- not allowed
  StepperIMEX_RK();

protected:

  std::string                                            description_;

  Teuchos::RCP<const RKButcherTableau<Scalar> >          explicitTableau_;
  Teuchos::RCP<const RKButcherTableau<Scalar> >          implicitTableau_;

  Scalar order_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >               stageX_;
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageF_;
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageG_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >               xTilde_;

  Teuchos::RCP<StepperObserver<Scalar> >         stepperObserver_;
  Teuchos::RCP<StepperIMEX_RKObserver<Scalar> >  stepperIMEX_RKObserver_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;

};


/** \brief Time-derivative interface for IMEX RK.
 *
 *  Given the stage state \f$X_i\f$ and
 *  \f[
 *    \tilde{X} = x_{n-1} +\Delta t \sum_{j=1}^{i-1} a_{ij}\,\dot{X}_{j},
 *  \f]
 *  compute the IMEX RK stage time-derivative,
 *  \f[
 *    \dot{X}_i = \frac{X_{i} - \tilde{X}}{a_{ii} \Delta t}
 *  \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperIMEX_RKTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperIMEX_RKTimeDerivative(
    Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  { initialize(s, xTilde); }

  /// Destructor
  virtual ~StepperIMEX_RKTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;

    // ith stage
    // s = 1/(dt*a_ii)
    // xOld = solution at beginning of time step
    // xTilde = xOld + dt*(Sum_{j=1}^{i-1} a_ij x_dot_j)
    // xDotTilde = - (s*x_i - s*xTilde)
    Thyra::V_StVpStV(xDot.ptr(),s_,*x,-s_,*xTilde_);
  }

  virtual void initialize(Scalar s,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  { s_ = s; xTilde_ = xTilde; }

private:

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde_;
  Scalar                                         s_;      // = 1/(dt*a_ii)
};


} // namespace Tempus
#endif // Tempus_StepperIMEX_RK_decl_hpp
