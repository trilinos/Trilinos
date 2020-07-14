// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExplicitRK_decl_hpp
#define Tempus_StepperExplicitRK_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperRKBase.hpp"
#include "Tempus_StepperExplicit.hpp"
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  #include "Tempus_StepperRKObserverComposite.hpp"
#endif


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
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Explicit RK with the application-action locations indicated.}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State $X \leftarrow x_{n-1}$ \Comment Set initial guess to last timestep.
 *    \For {$i = 0 \ldots s-1$}
 *        \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STAGE)}
 *        \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *      \If { i==0 and useFSAL and (previous step not failed) }
 *        \State tmp = $\dot{X}_0$
 *        \State $\dot{X}_0 = \dot{X}_s$
 *        \State $\dot{X}_s$ = tmp
 *        \State {\bf continue}
 *      \Else
 *        \State $X \leftarrow x_{n-1}
 *                  + \Delta t\,\sum_{j=1}^{i-1} a_{ij}\,\dot{X}_j$
 *        \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *        \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXPLICIT\_EVAL)}
 *        \State $\dot{X}_i \leftarrow \bar{f}(X_i,t_{n-1}+c_i\Delta t)$
 *      \EndIf
 *      \State {\it appAction.execute(solutionHistory, stepper, END\_STAGE)}
 *    \EndFor
 *    \State $x_n \leftarrow x_{n-1} + \Delta t\,\sum_{i=1}^{s}b_i\,\dot{X}_i$
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 *
 *   For Explicit RK, FSAL requires \f$c_1 = 0\f$, \f$c_s = 1\f$, and
 *   be stiffly accurate (\f$a_{sj} = b_j\f$).  An example of this is
 *   the Bogacki-Shampine 3(2) method.
 *   \f[
 *   \begin{array}{c|cccc}  0  & 0    &     &     &   \\
 *                         1/3 & 1/2  & 0   &     &   \\
 *                         2/3 & 0    & 3/4 & 0   &   \\
 *                          1  & 2/9  & 1/3 & 4/9 & 0 \\ \hline
 *                             & 2/9  & 1/3 & 4/9 & 0 \\
 *                             & 7/24 & 1/4 & 1/3 & 1/8 \end{array}
 *   \f]
 */
template<class Scalar>
class StepperExplicitRK : virtual public Tempus::StepperExplicit<Scalar>,
                          virtual public Tempus::StepperRKBase<Scalar>
{

public:

  /// \name Basic stepper methods
  //@{
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

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const;

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return false;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}

    void getValidParametersBasicERK(Teuchos::RCP<Teuchos::ParameterList> pl) const;
    virtual std::string getDescription() const = 0;
  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const;

  /// \name Accessors methods
  //@{
    /** \brief Use embedded if avialable. */
    virtual void setUseEmbedded(bool a) { useEmbedded_ = a; }
    virtual bool getUseEmbedded() const { return useEmbedded_; }
    virtual bool getUseEmbeddedDefault() const { return false; }
  //@}


protected:

  /// Default setup for constructor.
  virtual void setupDefault();

  /// Setup for constructor.
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  virtual void setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded);
#endif
  virtual void setup(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction);

  virtual void setupTableau() = 0;


  std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar> > > stageXDot_;

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  Teuchos::RCP<StepperRKObserverComposite<Scalar> >      stepperObserver_;
#endif

  // For Embedded RK
  bool useEmbedded_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               ee_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u0;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               abs_u;
  Teuchos::RCP<Thyra::VectorBase<Scalar> >               sc;

};

} // namespace Tempus

#endif // Tempus_StepperExplicitRK_decl_hpp
