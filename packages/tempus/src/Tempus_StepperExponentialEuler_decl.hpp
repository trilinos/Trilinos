// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperExponentialEuler_decl_hpp
#define Tempus_StepperExponentialEuler_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperExponentialEulerAppAction.hpp"


namespace Tempus {

/** \brief Exponential Euler time stepper.
 *
 *  For the explicit ODE system, \f$\dot{x} = \mathcal{F}(x,t)\f$,
 *  the solution, \f$x\f$, is determined using explicit evaluations of \f$\mathcal{F}\f$
 *  and exponentials of a linear operator \f$W\f$.
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Exponential Euler is:
 *
 *  \f{center}{
 *    \parbox{5in}{
 *    \rule{5in}{0.4pt} \\
 *    {\bf Algorithm} Backward Euler \\
 *    \rule{5in}{0.4pt} \vspace{-15pt}
 *    \begin{enumerate}
 *      \setlength{\itemsep}{0pt} \setlength{\parskip}{0pt} \setlength{\parsep}{0pt}
 *      \item {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *      \item {\bf Compute the explicit evaluation}  $f_{n-1} = \mathcal{F}(x_{n-1}, t_{n-1})$.
 *      \item {\it appAction.execute(solutionHistory, stepper, BEFORE\_EXP)}
 *      \item {\bf Compute $d_{n-1} = \exp(\Delta{t_{n-1}} W_{n-1}) f_{n-1}$.}
 *      \item {\it appAction.execute(solutionHistory, stepper, AFTER\_EXP)}
 *      \item $\dot{x}_n \leftarrow x_{n-1} + \Delta{t_n}d_n$
 *      \item {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *    \end{enumerate}
 *    \vspace{-10pt} \rule{5in}{0.4pt}
 *    }
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is not needed with Exponential Euler.
 *  The default is to set useFSAL=false, however useFSAL=true will also work
 *  but have no affect (i.e., no-op).
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \frac{\partial \mathcal{F}}{\partial x},
 *  \f]
 *  to obtain that from the implicit model evaluator, we set
 *  \f$ \alpha = 0 \f$ and \f$ \beta = 1 \f$.
 *  TODO: Need to use exponential model evaluator:
 *        explicit model evaluator plus W, or implicit model evaluator wrapper.
 */
template<class Scalar>
class StepperExponentialEuler :
    virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
  */
  StepperExponentialEuler();

  /// Constructor
  StepperExponentialEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> >& stepperEEAppAction);

  /// \name Basic stepper methods
  //@{
    virtual void setAppAction(
      Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > appAction);

    virtual Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > getAppAction() const
    { return stepperEEAppAction_; }

    /// Set the model
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel) override;

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState() override;
    virtual Scalar getOrder() const override {return 1.0;}
    virtual Scalar getOrderMin() const override {return 1.0;}
    virtual Scalar getOrderMax() const override {return 2.0;} //order is 2 for autonomous problems

    virtual bool isExplicit() const override {return false;}
    virtual bool isImplicit() const override {return true;}
    virtual bool isExplicitImplicit() const override
      {return isExplicit() && isImplicit();}
    virtual bool isOneStepMethod() const override {return true;}
    virtual bool isMultiStepMethod() const override {return !isOneStepMethod();}
    virtual OrderODE getOrderODE() const override {return FIRST_ORDER_ODE;}
  //@}

  // TODO: not sure what alpha should be for an exponential method
  // figure out where this is needed exterally (public) 
  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const override { return Scalar(1.0)/dt; }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar) const override { return Scalar(1.0); }

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const override;

private:

  /// Implementation of computeStep*() methods
  void computeStepResidDerivImpl(
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs,
    const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
    const Teuchos::Array<Scalar>& t,
    const Thyra::VectorBase<Scalar>& p,
    const int param_index,
    const int deriv_index = 0) const;

private:

  Teuchos::RCP<StepperExponentialEulerAppAction<Scalar> > stepperEEAppAction_;

};

/** \brief Time-derivative interface for Backward Euler.
 *
 *  Given the state \f$x\f$, compute the Backward Euler time-derivative,
 *  \f[
 *    \dot{x}_{n} = \frac{(x_{n} - x_{n-1})}{\Delta t_{n}}.
 *  \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperExponentialEulerTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:
  // TODO: remove this, not needed for exponential.
  /// Constructor
  StepperExponentialEulerTimeDerivative(
    Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld)
  { initialize(s, xOld); }

  /// Destructor
  virtual ~StepperExponentialEulerTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;
    // Calculate the Backward Euler x dot vector
    Thyra::V_StVpStV(xDot.ptr(),s_,*x,-s_,*xOld_);
  }

  virtual void initialize(Scalar s,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld)
  { s_ = s; xOld_ = xOld; }

private:

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld_;
  Scalar                                         s_;    // = 1.0/dt
};


/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template<class Scalar>
Teuchos::RCP<StepperExponentialEuler<Scalar> >
createStepperExponentialEuler(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl);


} // namespace Tempus

#endif // Tempus_StepperExponentialEuler_decl_hpp
