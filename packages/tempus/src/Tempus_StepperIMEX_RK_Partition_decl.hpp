// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperIMEX_RK_Partition_decl_hpp
#define Tempus_StepperIMEX_RK_Partition_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_Basic.hpp"
#include "Tempus_StepperIMEX_RKPartObserver.hpp"


namespace Tempus {

/** \brief Partitioned Implicit-Explicit Runge-Kutta (IMEX-RK) time stepper.
 *
 *  Partitioned IMEX-RK is similar to the IMEX-RK (StepperIMEX_RK),
 *  except a portion of the solution only requires explicit integration,
 *  and should not be part of the implicit solution to reduce computational
 *  costs.  Again our ODE can be written as
 *  \f{eqnarray*}{
 *    M(z,t)\, \dot{z} + G(z,t) + F(z,t) & = & 0, \\
 *    \mathcal{G}(\dot{z},z,t) + F(z,t) & = & 0,
 *  \f}
 *  but now
 *  \f[
 *    z     =\left\{\begin{array}{c} y\\ x \end{array}\right\},\;
 *    F(z,t)=\left\{\begin{array}{c} F^y(x,y,t)\\ F^x(x,y,t)\end{array}\right\},
 *    \mbox{ and }
 *    G(z,t)=\left\{\begin{array}{c} 0\\ G^x(x,y,t) \end{array}\right\}
 *  \f]
 *  where \f$z\f$ is the product vector of \f$y\f$ and \f$x\f$,
 *  \f$F(z,t)\f$ is still the "slow" physics (and evolved explicitly),
 *  and \f$G(z,t)\f$ is still the "fast" physics (and evolved implicitly),
 *  but a portion of the solution vector, \f$y\f$, is "explicit-only"
 *  and is only evolved by \f$F^y(x,y,t)\f$, while \f$x\f$ is the
 *  Implicit/Explicit (IMEX) solution vector, and is evolved explicitly by
 *  \f$F^x(x,y,t)\f$ evolved implicitly by \f$G^x(x,y,t)\f$.
 *  Note we can expand this to explicitly show all the terms as
 *  \f{eqnarray*}{
 *    & & M^y(x,y,t)\: \dot{y} + F^y(x,y,t) = 0, \\
 *    & & M^x(x,y,t)\: \dot{x} + F^x(x,y,t) + G^x(x,y,t) = 0, \\
 *  \f}
 *  or
 *  \f[
 *       \left\{ \begin{array}{c} \dot{y} \\ \dot{x} \end{array}\right\}
 *    +  \left\{ \begin{array}{c}    f^y  \\    f^x  \end{array}\right\}
 *    +  \left\{ \begin{array}{c}     0   \\    g^x  \end{array}\right\} = 0
 *  \f]
 *  where \f$f^y(x,y,t) = M^y(x,y,t)^{-1}\, F^y(x,y,t)\f$,
 *  \f$f^x(x,y,t) = M^x(x,y,t)^{-1}\, F^x(x,y,t)\f$, and
 *  \f$g^x(x,y,t) = M^x(x,y,t)^{-1}\, G^x(x,y,t)\f$,
 *  or
 *  \f[
 *    \dot{z} + f(x,y,t) + g(x,y,t) = 0,
 *  \f]
 *  where \f$f(x,y,t) = M(x,y,t)^{-1}\, F(x,y,t)\f$, and
 *  \f$g(x,y,t) = M(x,y,t)^{-1}\, G(x,y,t)\f$.
 *  Using Butcher tableaus for the
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
 *  the basic scheme for this partitioned, \f$s\f$-stage, IMEX-RK is
 *  \f[ \begin{array}{rcll}
 *   Z_i & = & Z_{n-1}
 *   - \Delta t \sum_{j=1}^{i-1} \hat{a}_{ij}\; f(Z_j,\hat{t}_j)
 *   - \Delta t \sum_{j=1}^i           a_{ij}\; g(Z_j,     t_j)
 *     &   \mbox{for } i=1\ldots s, \\
 *   z_n & = & z_{n-1}
 *   - \Delta t \sum_{i=1}^s \left[ \hat{b}_i\; f(Z_i,\hat{t}_i)
 *                                 +     b_i\;  g(Z_i,     t_i) \right] &
 *  \end{array} \f]
 *  or expanded
 *  \f[ \begin{array}{rcll}
 *   Y_i & = & y_{n-1}
 *   - \Delta t \sum_{j=1}^{i-1} \hat{a}_{ij}\; f^y(Z_j,\hat{t}_j)
 *      &   \mbox{for } i=1\ldots s,\\
 *   X_i & = & x_{n-1}
 *   - \Delta t \sum_{j=1}^{i-1} \hat{a}_{ij}\; f^x(Z_j,\hat{t}_j)
 *   - \Delta t \sum_{j=1}^i           a_{ij}\; g^x(Z_j,     t_j)
 *     &   \mbox{for } i=1\ldots s, \\
 *   y_n & = & y_{n-1}
 *   - \Delta t \sum_{i=1}^s \hat{b}_{i}\; f^y(X_i,Y_i,\hat{t}_i) & \\
 *   x_n & = & x_{n-1}
 *   - \Delta t \sum_{i=1}^s \left[ \hat{b}_i\; f^x(Z_i,\hat{t}_i)
 *                                 +     b_i\;  g^x(Z_i,     t_i) \right] &
 *  \end{array} \f]
 *  where \f$\hat{t}_i = t_{n-1}+\hat{c}_i\Delta t\f$ and
 *  \f$t_i = t_{n-1}+c_i\Delta t\f$.
 *
 *  For iterative solvers, it is useful to write the stage solutions as
 *  \f[
 *    Z_i = \tilde{Z} - a_{ii} \Delta t\, g(Z_i,t_i)
 *  \f]
 *  or expanded as
 *  \f[
 *    \left\{ \begin{array}{c}        Y_i \\        X_i  \end{array}\right\} =
 *    \left\{ \begin{array}{c} \tilde{Y}  \\ \tilde{X}_i \end{array}\right\}
 *    -  a_{ii} \Delta t
 *    \left\{ \begin{array}{c}        0   \\ g^x(Z_i,t_i) \end{array}\right\}
 *  \f]
 *  where
 *  \f{eqnarray*}{
 *    \tilde{Z} & = & z_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f(Z_j,\hat{t}_j) + a_{ij}\, g(Z_j, t_j)\right] \\
 *    \tilde{Y} & = & y_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f^y(Z_j,\hat{t}_j)\right] \\
 *    \tilde{X} & = & x_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f^x(Z_j,\hat{t}_j) +a_{ij}\, g^x(Z_j,t_j)\right] \\
 *  \f}
 *  and note that \f$Y_i = \tilde{Y}\f$.  Rearranging to solve for the
 *  implicit term
 *  \f{eqnarray*}{
 *    g  (Z_i,t_i) & = & - \frac{Z_i - \tilde{Z}}{a_{ii} \Delta t} \\
 *    g^x(Z_i,t_i) & = & - \frac{X_i - \tilde{X}}{a_{ii} \Delta t}
 *  \f}
 *  We additionally need the time derivative at each stage for the implicit
 *  solve.  Let us define the following time derivative for \f$x\f$ portion
 *  of the solution
 *  \f[
 *    \dot{X}_i(X_i,Y_i,t_i) + f^x(X_i,Y_i,t_i) + g^x(X_i,Y_i,t_i) = 0
 *  \f]
 *  where we split \f$Z_i\f$ arguments into \f$X_i\f$ and \f$Y_i\f$ to
 *  emphasize that \f$X_i\f$ is the solution for the implicit solve and
 *  \f$Y_i\f$ are parameters in this set of equations.  The above time
 *  derivative, \f$\dot{X}_i\f$, is NOT likely the same as the real time
 *  derivative, \f$\dot{x}(x(t_i), y(t_i), t_i)\f$, unless
 *  \f$\hat{c}_i = c_i \rightarrow \hat{t}_i = t_i\f$
 *  (Reasoning: \f$x(t_i) \neq X_i\f$ and \f$y(t_i) \neq Y_i\f$ unless
 *  \f$\hat{t}_i = t_i\f$).  Also note that the explicit term,
 *  \f$f^x(X_i,Y_i,t_i)\f$, is evaluated at the implicit stage time, \f$t_i\f$.
 *
 *  We can form the time derivative
 *  \f{eqnarray*}{
 *    \dot{X}(X_i,Y_i,t_i) & = & - g^x(X_i,Y_i,t_i) - f^x(X_i,Y_i,t_i) \\
 *    \dot{X}(X_i,Y_i,t_i) & = &
 *      \frac{X_i - \tilde{X}}{a_{ii} \Delta t} - f^x(X_i,Y_i,t_i) \\
 *  \f}
 *  Returning to the governing equation for the IMEX solution vector, \f$X_i\f$
 *  \f{eqnarray*}{
 *    M^x(X_i,Y_i,t_i)\,
 *    \dot{X}(X_i,Y_i,t_i) + F^x(X_i,Y_i,t_i) + G^x(X_i,Y_i,t_i) & = & 0 \\
 *    M^x(X_i,Y_i,t_i)\,
 *    \left[ \frac{X_i - \tilde{X}}{a_{ii} \Delta t} - f^x(X_i,Y_i,t_i) \right]
 *      + F^x(X_i,Y_i,t_i) + G^x(X_i,Y_i,t_i) & = & 0 \\
 *    M^x(X_i,Y_i,t_i)\,
 *    \left[ \frac{X_i - \tilde{X}}{a_{ii} \Delta t} \right]
 *    + G(X_i,Y_i,t_i) & = & 0 \\
 *  \f}
 *  Recall \f$\mathcal{G}^x(\dot{x},x,y,t) = M^x(x,y,t)\,\dot{x} + G^x(x,y,t)\f$
 *  and if we define a pseudo time derivative, which is equivalent to the
 *  time derivative for the implicit solve,
 *  \f[
 *    \tilde{\dot{X}} = \frac{X_i - \tilde{X}}{a_{ii} \Delta t},
 *  \f]
 *  we can write
 *  \f[
 *    \mathcal{G}^x(\tilde{\dot{X}},X_i,Y_i,t_i) =
 *       M^x(X_i,Y_i,t_i)\, \tilde{\dot{X}} + G^x(X_i,Y_i,t_i) = 0
 *  \f]
 *
 *  For general DIRK methods, we need to also handle the case when
 *  \f$a_{ii}=0\f$.  The IMEX stage values can be simply evaluated
 *  similiar to the "explicit-only" stage values, e.g.,
 *  \f[
 *     X_i = \tilde{X} = x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *            \hat{a}_{ij}\, f^x_j + a_{ij}\, g^x_j \right)
 *  \f]
 *  and then we can simply evaluate
 *  \f{eqnarray*}{
 *     f_i   & = & f  (Z_i,\hat{t}_i) \\
 *     g^x_i & = & g^x(X_i,Y_i,  t_i)
 *  \f}
 *  We can then form the time derivative as
 *  \f[
 *    \dot{X}_i = - g^x(X_i,Y_i,t_i) - f^x(X_i,Y_i,t_i)
 *  \f]
 *  but again note that the explicit term, \f$f^x(X_i,Y_i,t_i)\f$,
 *  is evaluated at the implicit stage time, \f$t_i\f$.
 *
 *  <b> Partitioned IMEX-RK Algorithm </b>
 *  The single-timestep algorithm for the partitioned IMEX-RK is
 *   - \f$Z_1 \leftarrow z_{n-1}\f$ (Recall \f$Z_i = \{Y_i,X_i\}^T\f$)
 *   - for \f$i = 1 \ldots s\f$ do
 *     - \f$Y_i = y_{n-1} -\Delta t \sum_{j=1}^{i-1} \hat{a}_{ij}\;f^y_j\f$
 *     - \f$\tilde{X} \leftarrow x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left[
 *            \hat{a}_{ij}\, f^x_j + a_{ij}\, g^x_j \right] \f$
 *     - if \f$a_{ii} = 0\f$
 *       - \f$X_i   \leftarrow \tilde{X}\f$
 *       - \f$g^x_i \leftarrow g^x(X_i,Y_i,t_i)\f$
 *     - else
 *       - Define \f$\tilde{\dot{X}}(X_i,Y_i,t_i)
 *            = \frac{X_i-\tilde{X}}{a_{ii} \Delta t}\f$
 *       - Solve \f$\mathcal{G}^x(\tilde{\dot{X}},X_i,Y_i,t_i) = 0\f$
 *         for \f$X_i\f$ where \f$Y_i\f$ are known parameters
 *       - \f$g^x_i \leftarrow - \tilde{\dot{X}}\f$
 *     - \f$f_i \leftarrow f(Z_i,\hat{t}_i)\f$
 *     - \f$\dot{X}_i
 *          = - g^x_i - f^x(Z_i,t_i)\f$ [Optional]
 *     - \f$\dot{Y}_i = - f^y_i \f$ [Optional]
 *   - end for
 *   - \f$z_n = z_{n-1} - \Delta t\,\sum_{i=1}^{s}\hat{b}_i\, f_i\f$
 *   - \f$x_n \mathrel{+{=}} - \Delta t\,\sum_{i=1}^{s} b_i\, g^x_i\f$
 *   - Solve \f$M(z_n) \dot{z}_n + F(z_n,t_n) + G(z_n,t_n) = 0\f$
 *       for \f$\dot{z}_n\f$ [Optional]
 *
 *  #### References
 *  -# Shadid, Cyr, Pawlowski, Widley, Scovazzi, Zeng, Phillips, Conde,
 *     Chuadhry, Hensinger, Fischer, Robinson, Rider, Niederhaus, Sanchez,
 *     "Towards an IMEX Monolithic ALE Method with Integrated UQ for
 *     Multiphysics Shock-hydro", SAND2016-11353, 2016, pp. 21-28.
 *  -# Cyr, "IMEX Lagrangian Methods", SAND2015-3745C.
 */
template<class Scalar>
class StepperIMEX_RK_Partition : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor to use default Stepper parameters.
  StepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    std::string stepperType = "Partitioned IMEX RK SSP2");

  /// Constructor to specialize Stepper parameters.
  StepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList);

  /// Constructor for StepperFactory.
  StepperIMEX_RK_Partition(
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
      const Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> > &
        modelPair);

    virtual void setModelPair(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel);

    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

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

  /// Pass initial guess to Newton solver (only relevant for implicit solvers)
  virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
     {initial_guess_ = initial_guess;}

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
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & Y,
    Scalar time, Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & G) const;

  void evalExplicitModel(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
    Scalar time, Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const;

private:

  /// Default Constructor -- not allowed
  StepperIMEX_RK_Partition();

protected:

  std::string                                            description_;
  Teuchos::RCP<const RKButcherTableau<Scalar> >          explicitTableau_;
  Teuchos::RCP<const RKButcherTableau<Scalar> >          implicitTableau_;

  Scalar order_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >               stageZ_;
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageF_;
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageGx_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >               xTilde_;

  Teuchos::RCP<StepperObserver<Scalar> >            stepperObserver_;
  Teuchos::RCP<StepperIMEX_RKPartObserver<Scalar> > stepperIMEX_RKPartObserver_;

  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;

};


/** \brief Time-derivative interface for Partitioned IMEX RK.
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
class StepperIMEX_RKPartTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperIMEX_RKPartTimeDerivative(
    Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  { initialize(s, xTilde); }

  /// Destructor
  virtual ~StepperIMEX_RKPartTimeDerivative() {}

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
#endif // Tempus_StepperIMEX_RK_Partition_decl_hpp
