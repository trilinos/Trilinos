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
#include "Tempus_StepperRKBase.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  #include "Tempus_StepperRKObserverComposite.hpp"
#endif


namespace Tempus {

/** \brief Implicit-Explicit Runge-Kutta (IMEX-RK) time stepper.
 *
 *  For the implicit ODE system, \f$ \mathcal{F}(\dot{x},x,t) = 0 \f$,
 *  we need to specialize this in order to separate the explicit, implicit,
 *  and temporal terms for the IMEX-RK time stepper,
 *  \f[
 *    M(x,t)\, \dot{x}(x,t) + G(x,t) + F(x,t) = 0
 *  \f]
 *  or
 *  \f[
 *    \mathcal{G}(\dot{x},x,t) + F(x,t) = 0,
 *  \f]
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
 *  explicit and implicit terms,
 *  \f[ \begin{array}{c|c}
 *    \hat{c} & \hat{a} \\ \hline
 *            & \hat{b}^T
 *  \end{array}
 *  \;\;\;\; \mbox{ and } \;\;\;\;
 *  \begin{array}{c|c}
 *    c & a \\ \hline
 *      & b^T
 *  \end{array}, \f]
 *  respectively, the basic IMEX-RK method for \f$s\f$-stages can be written as
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
 *  \f$t_i = t_{n-1}+c_i\Delta t\f$.  Note that the "slow" explicit physics,
 *  \f$f(X_j,\hat{t}_j)\f$, is evaluated at the explicit stage time,
 *  \f$\hat{t}_j\f$, and the "fast" implicit physics, \f$g(X_j,t_j)\f$, is
 *  evaluated at the implicit stage time, \f$t_j\f$.  We can write the stage
 *  solution, \f$X_i\f$, as
 *  \f[
 *    X_i = \tilde{X} - a_{ii} \Delta t\, g(X_i,t_i)
 *  \f]
 *  where
 *  \f[
 *    \tilde{X} = x_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *        \left(\hat{a}_{ij}\, f(X_j,\hat{t}_j) + a_{ij}\, g(X_j,t_j)\right)
 *  \f]
 *  Rewriting this in a form for Newton-type solvers, the implicit ODE is
 *  \f[
 *    \mathcal{G}(\tilde{\dot{X}},X_i,t_i) = \tilde{\dot{X}} + g(X_i,t_i) = 0
 *  \f]
 *  where we have defined a pseudo time derivative, \f$\tilde{\dot{X}}\f$,
 *  \f[
 *    \tilde{\dot{X}} \equiv \frac{X_i - \tilde{X}}{a_{ii} \Delta t}
 *    \quad \quad \left[ = -g(X_i,t_i)\right]
 *  \f]
 *  that can be used with the implicit solve but is <b>not</b> the stage
 *  time derivative, \f$\dot{X}_i\f$.  (Note that \f$\tilde{\dot{X}}\f$
 *  can be interpreted as the rate of change of the solution due to the
 *  implicit "fast" physics, and the "mass" version of the implicit ODE,
 *  \f$\mathcal{G}(\tilde{\dot{X}},X_i,t) = M(X_i,t_i)\, \tilde{\dot{X}}
 *  + G(X_i,t_i) = 0\f$, can also be used to solve for \f$\tilde{\dot{X}}\f$).
 *
 *  To obtain the stage time derivative, \f$\dot{X}_i\f$, we can evaluate
 *  the governing equation at the implicit stage time, \f$t_i\f$,
 *  \f[
 *    \dot{X}_i(X_i,t_i) + g(X_i,t_i) + f(X_i,t_i) = 0
 *  \f]
 *  Note that even the explicit term, \f$f(X_i,t_i)\f$, is evaluated at
 *  the implicit stage time, \f$t_i\f$.  Solving for \f$\dot{X}_i\f$, we find
 *  \f{eqnarray*}{
 *    \dot{X}(X_i,t_i) & = & - g(X_i,t_i) - f(X_i,t_i) \\
 *    \dot{X}(X_i,t_i) & = & \tilde{\dot{X}} - f(X_i,t_i)
 *  \f}
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n},
 *  \f]
 *  where \f$ \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}, \f$
 *  and \f$ \beta \equiv \frac{\partial x_n}{\partial x_n} = 1\f$. For the stage
 *  solutions, we are solving
 *  \f[
 *    \mathcal{G} = \tilde{\dot{X}} + g(X_i,t_i) = 0.
 *  \f]
 *  where \f$\mathcal{F}_n \rightarrow \mathcal{G}\f$,
 *  \f$x_n \rightarrow X_{i}\f$, and
 *  \f$\dot{x}_n(x_n) \rightarrow \tilde{\dot{X}}(X_{i})\f$.
 *  The time derivative for the implicit solves is
 *  \f[
 *    \tilde{\dot{X}} \equiv \frac{X_i - \tilde{X}}{a_{ii} \Delta t}
 *  \f]
 *  and we can determine that
 *  \f$ \alpha = \frac{1}{a_{ii} \Delta t} \f$
 *  and \f$ \beta = 1 \f$, and therefore write
 *  \f[
 *    W = \frac{1}{a_{ii} \Delta t}
 *        \frac{\partial \mathcal{G}}{\partial \tilde{\dot{X}}}
 *      + \frac{\partial \mathcal{G}}{\partial X_i}.
 *  \f]
 *
 *  <b>Explicit Stage in the Implicit Tableau.</b>
 *  For the special case of an explicit stage in the implicit tableau,
 *  \f$a_{ii}=0\f$, we find that the stage solution, \f$X_i\f$, is
 *  \f[
 *     X_i = x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *       \hat{a}_{ij}\,f(X_j,\hat{t}_j) + a_{ij}\,g(X_j,t_j) \right) = \tilde{X}
 *  \f]
 *  and the time derivative of the stage solution, \f$\dot{X}(X_i,t_i)\f$, is
 *  \f[
 *    \dot{X}_i(X_i,t_i) = - g(X_i,t_i) - f(X_i,t_i)
 *  \f]
 *  and again note that the explicit term, \f$f(X_i,t_i)\f$, is evaluated
 *  at the implicit stage time, \f$t_i\f$.
 *
 *  <b> IMEX-RK Algorithm </b>
 *
 *  The single-timestep algorithm for IMEX-RK is
 *
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{IMEX RK with the application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State $X \leftarrow x_{n-1}$ \Comment Set initial guess to last timestep.
 *    \For {$i = 0 \ldots s-1$}
 *      \State $\tilde{X} \leftarrow x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left(
 *            \hat{a}_{ij}\, f_j + a_{ij}\, g_j \right)$
 *      \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STAGE)}
 *      \State \Comment Implicit Tableau
 *      \If {$a_{ii} = 0$}
 *        \State $X \leftarrow \tilde{X}$
 *        \If {$a_{k,i} = 0 \;\forall k = (i+1,\ldots, s-1)$, $b(i) = 0$, $b^\ast(i) = 0$}
 *          \State $g_i \leftarrow 0$ \Comment{Not needed for later calculations.}
 *        \Else
 *          \State $g_i \leftarrow M(X, t_i)^{-1}\, G(X, t_i)$
 *        \EndIf
 *      \Else
 *        \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *        \If {``Zero initial guess.''}
 *          \State $X \leftarrow 0$
 *            \Comment{Else use previous stage value as initial guess.}
 *        \EndIf
 *        \State Solve $\mathcal{G}\left(\tilde{\dot{X}}
 *            = \frac{X-\tilde{X}}{a_{ii} \Delta t},X,t_i\right) = 0$ for $X$
 *        \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *        \State $\tilde{\dot{X}} \leftarrow \frac{X - \tilde{X}}{a_{ii} \Delta t}$
 *        \State $g_i \leftarrow - \tilde{\dot{X}}$
 *      \EndIf
 *      \State \Comment Explicit Tableau
 *      \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *      \State $f_i \leftarrow M(X,\hat{t}_i)^{-1}\, F(X,\hat{t}_i)$
 *      \State $\dot{X} \leftarrow - g_i - f_i$ [Optionally]
 *      \State {\it appAction.execute(solutionHistory, stepper, END\_STAGE)}
 *    \EndFor
 *    \State $x_n \leftarrow x_{n-1} - \Delta t\,\sum_{i=1}^{s}\hat{b}_i\,f_i
 *                                   - \Delta t\,\sum_{i=1}^{s}     b_i \,g_i$
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
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
 *  <tr><td> SSP1_111 <td> 1st
 *      <td> \f[ \begin{array}{c|c}
 *             1 & 1 \\ \hline
 *               & 1
 *           \end{array} \f]
 *      <td> \f[ \begin{array}{c|c}
 *             0 & 0 \\ \hline
 *               & 1
 *           \end{array} \f]
 *  <tr><td> IMEX RK SSP2\n SSP2_222_L\n \f$\gamma = 1-1/\sqrt{2}\f$ <td> 2nd
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
 *  <tr><td> SSP2_222\n SSP2_222_A\n \f$\gamma = 1/2\f$ <td> 2nd
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
 *  <tr><td> IMEX RK SSP3\n SSP3_332\n \f$\gamma = 1/ (2 + \sqrt{2})\f$ <td> 3rd
 *      <td> \f[ \begin{array}{c|ccc}
 *             0  &  0  &     &     \\
 *             1  &  1  &  0  &     \\
 *            1/2 & 1/4 & 1/4 &  0  \\ \hline
 *                & 1/6 & 1/6 & 4/6  \end{array} \f]
 *      <td> \f[ \begin{array}{c|ccc}
 *              \gamma  & \gamma      &         &        \\
 *             1-\gamma & 1-2\gamma   & \gamma  &        \\
 *             1-2      & 1/2 -\gamma & 0       & \gamma \\ \hline
 *                      & 1/6         & 1/6     & 2/3    \end{array} \f]
 *  <tr><td> IMEX RK ARS 233\n ARS 233\n \f$\gamma = (3+\sqrt{3})/6\f$ <td> 3rd
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
 *  The First-Step-As-Last (FSAL) principle is not valid for IMEX RK.
 *  The default is to set useFSAL=false, and useFSAL=true will result
 *  in an error.
 *
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
 *
 */
template<class Scalar>
class StepperIMEX_RK : virtual public Tempus::StepperImplicit<Scalar>,
                       virtual public Tempus::StepperRKBase<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
  */
  StepperIMEX_RK();

  /// Constructor for all member data.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperObserver<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    std::string stepperType,
    Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau,
    Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau,
    Scalar order);
#endif
  StepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    std::string stepperType,
    Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau,
    Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau,
    Scalar order);


  /// \name Basic stepper methods
  //@{
    /// Returns the explicit tableau!
    virtual Teuchos::RCP<const RKButcherTableau<Scalar> > getTableau() const
    { return getExplicitTableau(); }

    /// Set both the explicit and implicit tableau from ParameterList
    virtual void setTableaus( std::string stepperType = "",
      Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau = Teuchos::null,
      Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau = Teuchos::null);

    /// Return explicit tableau.
    virtual Teuchos::RCP<const RKButcherTableau<Scalar> > getExplicitTableau() const
    { return explicitTableau_; }

    /// Set the explicit tableau from tableau
    virtual void setExplicitTableau(
      Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau);

    /// Return implicit tableau.
    virtual Teuchos::RCP<const RKButcherTableau<Scalar> > getImplicitTableau() const
    { return implicitTableau_; }

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

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual Teuchos::RCP<StepperObserver<Scalar> > getObserver() const
    { return this->stepperObserver_; }
#endif
    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

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

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
  //@}

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > >& getStageF() {return stageF_;};
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > >& getStageG() {return stageG_;};
  Teuchos::RCP<Thyra::VectorBase<Scalar> >& getXTilde() {return xTilde_;};

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const
  {
    const Teuchos::SerialDenseMatrix<int,Scalar> & A = implicitTableau_->A();
    return Scalar(1.0)/(dt*A(0,0));  // Getting the first diagonal coeff!
  }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar   ) const { return Scalar(1.0); }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const;

  void evalImplicitModelExplicitly(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
    Scalar time, Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & G) const;

  void evalExplicitModel(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X,
    Scalar time, Scalar stepSize, Scalar stageNumber,
    const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const;

  virtual bool getICConsistencyCheckDefault() const { return false; }

  void setOrder(Scalar order) { order_ = order; }

protected:

  Teuchos::RCP<const RKButcherTableau<Scalar> >          explicitTableau_;
  Teuchos::RCP<const RKButcherTableau<Scalar> >          implicitTableau_;

  Scalar order_;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageF_;
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageG_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >               xTilde_;

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  Teuchos::RCP<StepperRKObserverComposite<Scalar> >        stepperObserver_;
#endif

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
