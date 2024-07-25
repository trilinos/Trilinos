//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperExplicitRK_decl_hpp
#define Tempus_StepperExplicitRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKBase.hpp"
#include "Tempus_StepperExplicit.hpp"

namespace Tempus {

/** \brief Explicit Runge-Kutta time stepper.
 *
 *  For the explicit ODE system,
 *  \f[
 *    \dot{x} = \bar{f}(x,t),
 *  \f]
 *  the general explicit Runge-Kutta method for \f$s\f$-stages can be
 *  written as
 *  \f[
 *    X_{i} = x_{n-1}
 *    + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\bar{f}(X_{j},t_{n-1}+c_{j}\Delta t)
 *  \f]
 *  \f[
 *    x_{n} = x_{n-1}
 *    + \Delta t\,\sum_{i=1}^{s}b_{i}\,\bar{f}(X_{i},t_{n-1}+c_{i}\Delta t)
 *  \f]
 *  where \f$X_{i}\f$ are intermediate approximations to the solution
 *  at times, \f$t_{n-1}+c_{i}\Delta t\f$, (stage solutions) which may
 *  be correct to a lower order of accuracy than the solution, \f$x_{n}\f$.
 *  We should note that these lower-order approximations are combined
 *  through \f$b_{i}\f$ so that error terms cancel out and produce a
 *  more accurate solution. Note for explicit RK that \f$a_{ij}=0\f$ for
 *  \f$j \leq i\f$ and does not require any solves.
 *  Note that the stage time derivatives are
 *  \f[
 *    \dot{X}_{i} = \bar{f}(X_{i},t_{n-1}+c_{i}\Delta t),
 *  \f]
 *  and the time derivative by definition is
 *  \f[
 *    \dot{x}_{n} = \bar{f}(x_{n},t_{n}),
 *  \f]
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Explicit RK is
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Explicit RK \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt}
 * \setlength{\parsep}{0pt} \item $X \leftarrow x_{n-1}$ \item {\it
 * appAction.execute(solutionHistory, stepper, BEGIN\_STEP)} \item {\bf for {$i
 * = 0 \ldots s-1$}} \item \quad $X \leftarrow x_{n-1}
 *                    + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_j$
 *      \item \quad {\it appAction.execute(solutionHistory, stepper,
 * BEGIN\_STAGE)} \item \quad {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_SOLVE)} \item \quad {\it appAction.execute(solutionHistory, stepper,
 * AFTER\_SOLVE)} \item \quad {\it appAction.execute(solutionHistory, stepper,
 * BEFORE\_EXPLICIT\_EVAL)} \item \quad {\bf if (i=0 and useFSAL and (previous
 * step not failed)) then} \item \qquad  tmp = $\dot{X}_0$ \item \qquad
 * $\dot{X}_0 = \dot{X}_s$ \item \qquad  $\dot{X}_s$ = tmp \item \qquad  {\bf
 * continue} \item \quad {\bf else} \item \qquad  $\dot{X}_i \leftarrow
 * \bar{f}(X_i,t_{n-1}+c_i\Delta t)$ \item \quad {\bf endif} \item \quad {\it
 * appAction.execute(solutionHistory, stepper, END\_STAGE)} \item {\bf end for}
 *      \item $x_n \leftarrow x_{n-1} + \Delta t\,\sum_{i=1}^{s}b_i\,\dot{X}_i$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *   For Explicit RK, FSAL requires \f$c_1 = 0\f$, \f$c_s = 1\f$, and
 *   be stiffly accurate (\f$a_{sj} = b_j\f$).  An example of this is
 *   the Bogacki-Shampine 3(2) method.
 *   \f[
 *   \begin{array}{c|cccc}  0  & 0    &     &     &   \\
 *                         1/2 & 1/2  & 0   &     &   \\
 *                         3/4 & 0    & 3/4 & 0   &   \\
 *                          1  & 2/9  & 1/3 & 4/9 & 0 \\ \hline
 *                             & 2/9  & 1/3 & 4/9 & 0 \\
 *                             & 7/24 & 1/4 & 1/3 & 1/8 \end{array}
 *   \f]
 */
template <class Scalar>
class StepperExplicitRK : virtual public Tempus::StepperExplicit<Scalar>,
                          virtual public Tempus::StepperRKBase<Scalar> {
 public:
  /// \name Basic stepper methods
  //@{
  /// Initialize during construction and after changing input parameters.
  virtual void initialize();

  /// Set model
  virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);

  /// Set the initial conditions and make them consistent.
  virtual void setInitialConditions(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Take the specified timestep, dt, and return true if successful.
  virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  /// Get a default (initial) StepperState
  virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  virtual Scalar getInitTimeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const;

  virtual bool isExplicit() const { return true; }
  virtual bool isImplicit() const { return false; }
  virtual bool isExplicitImplicit() const
  {
    return isExplicit() && isImplicit();
  }
  virtual bool isOneStepMethod() const { return true; }
  virtual bool isMultiStepMethod() const { return !isOneStepMethod(); }

  virtual OrderODE getOrderODE() const { return FIRST_ORDER_ODE; }

  virtual std::string getDescription() const = 0;
  //@}

  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<Teuchos::ParameterList> getValidParametersBasicERK() const;

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream& out) const;

 protected:
  /// Default setup for constructor.
  virtual void setupDefault();

  /// Setup for constructor.
  virtual void setup(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
      bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
      bool useEmbedded,
      const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction);

  virtual void setupTableau() = 0;

  virtual void setEmbeddedMemory();

  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_;
};

}  // namespace Tempus

#endif  // Tempus_StepperExplicitRK_decl_hpp
