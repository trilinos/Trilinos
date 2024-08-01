//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperIMEX_RK_Partition_decl_hpp
#define Tempus_StepperIMEX_RK_Partition_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKBase.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluatorPairPartIMEX_Basic.hpp"

namespace Tempus {

/** \brief Partitioned Implicit-Explicit Runge-Kutta (IMEX-RK) time stepper.
 *
 *  Partitioned IMEX-RK is similar to the IMEX-RK (StepperIMEX_RK),
 *  except a portion of the solution only requires explicit integration,
 *  and should not be part of the implicit solution to reduce computational
 *  costs.  Again our implicit ODE can be written as
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
 *  \f$F^x(x,y,t)\f$ and is evolved implicitly by \f$G^x(x,y,t)\f$.
 *  Note we can expand this to show all the terms as
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
 *  where
 *  \f{eqnarray*}{
 *    f^y(x,y,t) & = & M^y(x,y,t)^{-1}\, F^y(x,y,t), \\
 *    f^x(x,y,t) & = & M^x(x,y,t)^{-1}\, F^x(x,y,t), \\
 *    g^x(x,y,t) & = & M^x(x,y,t)^{-1}\, G^x(x,y,t),
 *  \f}
 *  or
 *  \f[
 *    \dot{z} + g(x,y,t) + f(x,y,t) = 0,
 *  \f]
 *  where \f$f(x,y,t) = M(x,y,t)^{-1}\, F(x,y,t)\f$, and
 *  \f$g(x,y,t) = M(x,y,t)^{-1}\, G(x,y,t)\f$.
 *  Using Butcher tableaus for the explicit and implicit terms
 *  \f[ \begin{array}{c|c}
 *    \hat{c} & \hat{a} \\ \hline
 *            & \hat{b}^T
 *  \end{array}
 *  \;\;\;\; \mbox{ and } \;\;\;\;
 *  \begin{array}{c|c}
 *    c & a \\ \hline
 *      & b^T
 *  \end{array}, \f]
 *  respectively, the basic scheme for this partitioned, \f$s\f$-stage,
 *  IMEX-RK method is
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
 *   - \Delta t \sum_{i=1}^s \hat{b}_{i}\; f^y(Z_i,\hat{t}_i) & \\
 *   x_n & = & x_{n-1}
 *   - \Delta t \sum_{i=1}^s \left[ \hat{b}_i\; f^x(Z_i,\hat{t}_i)
 *                                 +     b_i\;  g^x(Z_i,     t_i) \right] &
 *  \end{array} \f]
 *  where \f$\hat{t}_i = t_{n-1}+\hat{c}_i\Delta t\f$ and
 *  \f$t_i = t_{n-1}+c_i\Delta t\f$.  Note that the "slow" explicit physics,
 *  \f$f^y(Z_j,\hat{t}_j)\f$ and \f$f^x(Z_j,\hat{t}_j)\f$, is evaluated at
 *  the explicit stage time, \f$\hat{t}_j\f$, and the "fast" implicit physics,
 *  \f$g^x(Z_j,t_j)\f$, is evaluated at the implicit stage time, \f$t_j\f$.
 *  We can write the stage solution, \f$Z_i\f$, as
 *  \f[
 *    Z_i = \tilde{Z} - a_{ii} \Delta t\, g(Z_i,t_i)
 *  \f]
 *  where
 *  \f[
 *    \tilde{Z} = z_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f(Z_j,\hat{t}_j) + a_{ij}\, g(Z_j, t_j)\right]
 *  \f]
 *  or in expanded form as
 *  \f[
 *    \left\{ \begin{array}{c}        Y_i \\        X_i  \end{array}\right\} =
 *    \left\{ \begin{array}{c} \tilde{Y}  \\ \tilde{X}_i \end{array}\right\}
 *    -  a_{ii} \Delta t
 *    \left\{ \begin{array}{c}        0   \\ g^x(Z_i,t_i) \end{array}\right\}
 *  \f]
 *  where
 *  \f{eqnarray*}{
 *    \tilde{Y} & = & y_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f^y(Z_j,\hat{t}_j)\right] \\
 *    \tilde{X} & = & x_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f^x(Z_j,\hat{t}_j) +a_{ij}\, g^x(Z_j,t_j)\right] \\
 *  \f}
 *  and note that \f$Y_i = \tilde{Y}\f$.
 *
 *  Noting that we will be solving the implicit ODE for \f$\dot{X}_i\f$,
 *  we can write
 *  \f[
 *    \mathcal{G}^x(\tilde{\dot{X}},X_i,Y_i,t_i) =
 *      \tilde{\dot{X}} + g^x(X_i,Y_i,t_i) = 0
 *  \f]
 *  where we have defined a pseudo time derivative, \f$\tilde{\dot{X}}\f$,
 *  \f[
 *    \tilde{\dot{X}} \equiv \frac{X_i - \tilde{X}}{a_{ii} \Delta t}
 *    \quad \quad \left[ = -g^x(X_i,Y_i,t_i)\right]
 *  \f]
 *  that can be used with the implicit solve but is <b>not</b> the stage
 *  time derivative for the IMEX equations, \f$\dot{X}_i\f$.
 *  (Note that \f$\tilde{\dot{X}}\f$
 *  can be interpreted as the rate of change of the solution due to the
 *  implicit "fast" physics.)
 *  Note that we are solving for \f$X_i\f$, and \f$Y_i\f$ are included as
 *  parameters possibly needed in the IMEX equations.
 *
 *  To obtain the stage time derivative, \f$\dot{Z}_i\f$, for the
 *  entire system, we can evaluate
 *  the governing equation at the implicit stage time, \f$t_i\f$,
 *  \f[
 *    \dot{Z}_i(Z_i,t_i) + f(Z_i,t_i) + g(Z_i,t_i) = 0
 *  \f]
 *  The above time derivative, \f$\dot{Z}_i\f$, is likely <b>not</b>
 *  the same as the real time derivative, \f$\dot{x}(x(t_i), y(t_i), t_i)\f$,
 *  unless \f$\hat{c}_i = c_i \rightarrow \hat{t}_i = t_i\f$
 *  (Reasoning: \f$x(t_i) \neq X_i\f$ and \f$y(t_i) \neq Y_i\f$ unless
 *  \f$\hat{t}_i = t_i\f$).  Also note that the explicit term,
 *  \f$f(Z_i,t_i)\f$, is evaluated at the implicit stage time, \f$t_i\f$.
 *  Solving for \f$\dot{Z}_i\f$, we find
 *  \f[
 *    \dot{Z}(Z_i,t_i) = - g(Z_i,t_i) - f(Z_i,t_i)
 *  \f]
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n},
 *  \f]
 *  where \f$ \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}, \f$
 *  and \f$ \beta \equiv \frac{\partial x_n}{\partial x_n} = 1\f$. For the
 *  IMEX equations, we are solving
 *  \f[
 *    \mathcal{G}^x(\tilde{\dot{X}},X_i,Y_i,t_i) =
 *      \tilde{\dot{X}} + g^x(X_i,Y_i,t_i) = 0
 *  \f]
 *  where \f$\mathcal{F}_n \rightarrow \mathcal{G}^x\f$,
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
 *        \frac{\partial \mathcal{G}^x}{\partial \tilde{\dot{X}}}
 *      + \frac{\partial \mathcal{G}^x}{\partial X_i}.
 *  \f]
 *
 *  <b>Explicit Stage in the Implicit Tableau.</b>
 *  For general DIRK methods, we need to also handle the case when
 *  \f$a_{ii}=0\f$.  The IMEX stage values can be simply evaluated
 *  similiar to the "explicit-only" stage values, e.g.,
 *  \f[
 *     X_i = \tilde{X} = x_{n-1} - \Delta t \sum_{j=1}^{i-1}
 *      \left[\hat{a}_{ij}\, f^x(Z_j,\hat{t}_j) +a_{ij}\, g^x(Z_j,t_j)\right]
 *  \f]
 *  and the time derivative of the stage solution is
 *  \f[
 *    \dot{X}_i = - g^x(X_i,Y_i,t_i) - f^x(X_i,Y_i,t_i)
 *  \f]
 *  but again note that the explicit term, \f$f^x(X_i,Y_i,t_i)\f$,
 *  is evaluated at the implicit stage time, \f$t_i\f$.
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for the partitioned IMEX-RK is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Partitioned IMEX-RK \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item $Z \leftarrow z_{n-1}$ \hfill {\it * Recall
 * $Z_i = \{Y_i,X_i\}^T$} \item {\it appAction.execute(solutionHistory, stepper,
 * BEGIN\_STEP)} \item {\bf for ($i = 0 \ldots s-1$)} \item \quad  $Y_i =
 * y_{n-1} -\Delta t \sum_{j=1}^{i-1} \hat{a}_{ij}\;f^y_j$ \item \quad
 * $\tilde{X} \leftarrow x_{n-1} - \Delta t\,\sum_{j=1}^{i-1} \left[
 *                                       \hat{a}_{ij}\, f^x_j + a_{ij}\, g^x_j
 * \right]$ \item \quad  {\it appAction.execute(solutionHistory, stepper,
 * BEGIN\_STAGE)} \item \quad  \hfill {\bf Implicit Tableau} \item \quad  {\bf
 * if ($a_{ii} = 0$) then} \item \qquad    $X_i   \leftarrow \tilde{X}$ \item
 * \qquad    {\bf if ($a_{k,i} = 0 \;\forall k = (i+1,\ldots, s-1)$, $b(i) = 0$,
 * $b^\ast(i) = 0$) then} \item \qquad \quad  $g^x_i \leftarrow 0$ \hfill {\it *
 * Not needed for later calculations.} \item \qquad    {\bf else} \item \qquad
 * \quad  $g^x_i \leftarrow g^x(X_i,Y_i,t_i)$ \item \qquad    {\bf endif} \item
 * \quad  {\bf else} \item \qquad  {\it appAction.execute(solutionHistory,
 * stepper, BEFORE\_SOLVE)} \item \qquad  {\bf if (``Zero initial guess.'')
 * then} \item \qquad \quad   $X \leftarrow 0$ \hfill {\it * Else use previous
 * stage value as initial guess.} \item \qquad  {\bf endif} \item \qquad  {\bf
 * Solve $\mathcal{G}^x\left(\tilde{\dot{X}} = \frac{X-\tilde{X}}{a_{ii} \Delta
 * t},X,Y,t_i\right) = 0$ for $X$ where $Y$ is known.} \item \qquad   {\it
 * appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)} \item \qquad
 * $\tilde{\dot{X}} \leftarrow \frac{X - \tilde{X}}{a_{ii} \Delta t}$ \item
 * \qquad   $g_i \leftarrow - \tilde{\dot{X}}$ \item \quad  {\bf endif} \item
 * \quad  \hfill {\bf Explicit Tableau} \item \quad  {\it
 * appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)} \item
 * \quad  $f_i \leftarrow M(X,\hat{t}_i)^{-1}\, F(X,\hat{t}_i)$ \item \quad
 * $\dot{Z} \leftarrow - g(Z_i,t_i) - f(Z_i,t_i)$ [Optionally] \item \quad  {\it
 * appAction.execute(solutionHistory, stepper, END\_STAGE)} \item {\bf end for}
 *      \item $z_n = z_{n-1} - \Delta t\,\sum_{i=1}^{s}\hat{b}_i\, f_i$
 *      \item $x_n \mathrel{+{=}} - \Delta t\,\sum_{i=1}^{s} b_i\, g^x_i$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The following table contains the pre-coded IMEX-RK tableaus.
 *  <table>
 *  <caption id="multi_row_part">Partitioned IMEX-RK Tableaus</caption>
 *  <tr><th> Name  <th> Order <th> Implicit Tableau <th> Explicit Tableau
 *  <tr><td> Partitioned IMEX RK 1st order  <td> 1st
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
 *  <tr><td> Partitioned IMEX RK SSP2\n \f$\gamma = 1-1/\sqrt{2}\f$ <td> 2nd
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
 *  <tr><td> Partitioned IMEX RK ARS 233\n \f$\gamma = (3+\sqrt{3})/6\f$ <td>
 * 3rd <td> \f[ \begin{array}{c|ccc}
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
 *  The First-Same-As-Last (FSAL) principle is not valid for IMEX RK Partition.
 *  The default is to set useFSAL=false, and useFSAL=true will result
 *  in a warning.
 *
 *  #### References
 *  -# Shadid, Cyr, Pawlowski, Widley, Scovazzi, Zeng, Phillips, Conde,
 *     Chuadhry, Hensinger, Fischer, Robinson, Rider, Niederhaus, Sanchez,
 *     "Towards an IMEX Monolithic ALE Method with Integrated UQ for
 *     Multiphysics Shock-hydro", SAND2016-11353, 2016, pp. 21-28.
 *  -# Cyr, "IMEX Lagrangian Methods", SAND2015-3745C.
 */
template <class Scalar>
class StepperIMEX_RK_Partition : virtual public Tempus::StepperImplicit<Scalar>,
                                 virtual public Tempus::StepperRKBase<Scalar> {
 public:
  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
   */
  StepperIMEX_RK_Partition(
      std::string stepperType = "Partitioned IMEX RK SSP2");

  /// Constructor to for all member data.
  StepperIMEX_RK_Partition(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
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
  {
    return getExplicitTableau();
  }

  /// Set both the explicit and implicit tableau from ParameterList
  virtual void setTableaus(
      std::string stepperType = "",
      Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau =
          Teuchos::null,
      Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau =
          Teuchos::null);

  virtual void setTableausPartition(Teuchos::RCP<Teuchos::ParameterList> pl,
                                    std::string stepperType);

  /// Return explicit tableau.
  virtual Teuchos::RCP<const RKButcherTableau<Scalar> > getExplicitTableau()
      const
  {
    return explicitTableau_;
  }

  /// Set the explicit tableau from tableau
  virtual void setExplicitTableau(
      Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau);

  /// Return implicit tableau.
  virtual Teuchos::RCP<const RKButcherTableau<Scalar> > getImplicitTableau()
      const
  {
    return implicitTableau_;
  }

  /// Set the implicit tableau from tableau
  virtual void setImplicitTableau(
      Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau);

  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const
  {
    return this->wrapperModel_;
  }

  virtual void setModelPair(
      const Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >&
          modelPair);

  virtual void setModelPair(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& explicitModel,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& implicitModel);

  /// Initialize during construction and after changing input parameters.
  virtual void initialize();

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getOrder() const { return order_; }
  virtual Scalar getOrderMin() const { return order_; }
  virtual Scalar getOrderMax() const { return order_; }

  virtual bool isExplicit() const { return true; }
  virtual bool isImplicit() const { return true; }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const { return true; }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }
  virtual OrderODE getOrderODE() const { return FIRST_ORDER_ODE; }
  //@}

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > >& getStageF()
  {
    return stageF_;
  }
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > >& getStageGx()
  {
    return stageGx_;
  }
  Teuchos::RCP<Thyra::VectorBase<Scalar> >& getXTilde() { return xTilde_; }

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const
  {
    const Teuchos::SerialDenseMatrix<int, Scalar>& A = implicitTableau_->A();
    return Scalar(1.0) / (dt * A(0, 0));  // Getting the first diagonal coeff!
  }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta(const Scalar) const { return Scalar(1.0); }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

  void evalImplicitModelExplicitly(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& X,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& Y, Scalar time,
      Scalar stepSize, Scalar stageNumber,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& G) const;

  void evalExplicitModel(
      const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& X, Scalar time,
      Scalar stepSize, Scalar stageNumber,
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >& F) const;

  void setOrder(Scalar order) { order_ = order; }

 protected:
  Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau_;
  Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau_;

  Scalar order_;

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageF_;
  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageGx_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> > xTilde_;
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
  : virtual public Tempus::TimeDerivative<Scalar> {
 public:
  /// Constructor
  StepperIMEX_RKPartTimeDerivative(
      Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  {
    initialize(s, xTilde);
  }

  /// Destructor
  virtual ~StepperIMEX_RKPartTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
      Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDot,
      Teuchos::RCP<Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;

    // ith stage
    // s = 1/(dt*a_ii)
    // xOld = solution at beginning of time step
    // xTilde = xOld + dt*(Sum_{j=1}^{i-1} a_ij x_dot_j)
    // xDotTilde = - (s*x_i - s*xTilde)
    Thyra::V_StVpStV(xDot.ptr(), s_, *x, -s_, *xTilde_);
  }

  virtual void initialize(Scalar s,
                          Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde)
  {
    s_      = s;
    xTilde_ = xTilde;
  }

 private:
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xTilde_;
  Scalar s_;  // = 1/(dt*a_ii)
};

/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperIMEX_RK_Partition<Scalar> > createStepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType, Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus
#endif  // Tempus_StepperIMEX_RK_Partition_decl_hpp
