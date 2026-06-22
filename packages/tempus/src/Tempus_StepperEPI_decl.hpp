// @HEADER
// ****************************************************************************
// TODO
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperEPI_decl_hpp
#define Tempus_StepperEPI_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperEPIAppAction.hpp"

#include "Tempus_PhiEvaluator.hpp"


template <class Scalar>
class ExponentialODEParameters {
 public:
  /// Constructor
  ExponentialODEParameters() : timeStepSize_(Scalar(0.0)), stageNumber_(0) {}

  /// Constructor
  ExponentialODEParameters(Scalar timeStepSize, int stageNumber = 0)
    : timeStepSize_(timeStepSize), stageNumber_(stageNumber)
  {
  }

  Scalar timeStepSize_;
  int stageNumber_;
};


namespace Tempus {
// TODO: FURKAN: FIX THE DESCRIPTION
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
 *    {\bf Algorithm} Exponential Euler \\
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
class StepperEPI :
    virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
  */
  StepperEPI();

  /// Constructor
  StepperEPI(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction);

  /// \name Basic stepper methods
  //@{
    virtual void setAppAction(
      Teuchos::RCP<StepperEPIAppAction<Scalar> > appAction);

    virtual Teuchos::RCP<StepperEPIAppAction<Scalar> > getAppAction() const
    { return stepperEPIAppAction_; }

    /// Set the model
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel) override;

    /// Set the order
    void setOrder(Scalar order) {this->order_ = order;}

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions(
        const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    /// Set the initial conditions and make them consistent.
    virtual void evaluateExponentialODE(Teuchos::RCP<Thyra::VectorBase<Scalar> >& f,
        const Teuchos::RCP<Thyra::VectorBase<Scalar> >& x,
        const Teuchos::RCP<Thyra::VectorBase<Scalar> >& xDot, const Scalar time,
        const Teuchos::RCP<ExponentialODEParameters<Scalar> >& p);

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) override;

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState() override;
    virtual Scalar getOrder() const override {return order_;}
    virtual Scalar getOrderMin() const override {return 2.0;}
    virtual Scalar getOrderMax() const override {return 3.0;}
    virtual void setUseFSAL(bool a) override
    {
      this->useFSAL_       = a;
      this->isInitialized_ = false;
    }

    virtual bool isExplicit() const override {return false;}

    /// Get the implicit/explicit type: we return true, since we rely on the implicit ModelEvaluator
    virtual bool isImplicit() const override {return true;}
    virtual bool isExplicitImplicit() const override
      {return isExplicit() && isImplicit();}
    virtual bool isOneStepMethod() const override {return false;}
    virtual bool isMultiStepMethod() const override {return !isOneStepMethod();}
    virtual OrderODE getOrderODE() const override {return FIRST_ORDER_ODE;}
  //@}

  // TODO: not sure what alpha should be for an exponential method
  // figure out where this is needed exterally (public)
  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const override { return Scalar(1.0)/dt; }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar) const override { return Scalar(1.0); }

  /// Set StepperExponential member data from the ParameterList.
  void setStepperExponentialValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const override;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const override;

  void setPhiEvaluatorParameterList(const Teuchos::RCP<Teuchos::ParameterList>& pl)
  { phiEvaluatorPL_ = pl; }

 private:
  /// Compute the temporal finite difference dt_Mf_deriv
  ///   d/dt (-M * F(x,t))
  void computeTemporalFD(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf
  );

  /// Compute the nonlinear remainder:
  ///   remf = -M * (F(xr,tr) - F(x0,t0) - J_{x0} * (xr-x0) - F'(t0) * (tr-t0))
  /// including multiple of negative mass matrix (-M).
  ///
  /// dt is the current time-step, not necessarily (tr-t0)
  /// Mf contains already evaluated -M*F(x0,t0)
  /// dt_Mf_deriv contains already evaluated dt*M*F'(t0)
  /// Mfr is unused, currently, can contain -M*F(xr,tr)
  void computeRemf(
    Teuchos::RCP<Thyra::VectorBase<Scalar>>& remf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& xr,
    const Scalar tr,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x0,
    const Scalar t0,
    const Scalar dt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mf,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& dt_Mf_deriv,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& Mfr = Teuchos::null
  );

private:

  Teuchos::RCP<StepperEPIAppAction<Scalar> > stepperEPIAppAction_;

  Teuchos::RCP<PhiEvaluator<Scalar> > phiEvaluator_;

  Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL_;

  /// Finite difference step size used for RHS time derivative estimation
  /// needed for nonautonomous correction.
  Scalar temporal_finite_difference_eps_;

  /// temporal Integration order
  Scalar order_;

  /// Number of time steps to wait between adapt PhiEvaluator calls
  int adapt_phi_evaluator_interval_;
};


/// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template<class Scalar>
Teuchos::RCP<StepperEPI<Scalar> >
createStepperEPI(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  Teuchos::RCP<Teuchos::ParameterList> pl);


// ----------------------------------------------------------------------------
/** \brief EPI2 Definition.
 *
 *  See Tempus_StepperEPI for additional details.
 */
template <class Scalar>
class StepperExponential_EPI2 : virtual public StepperEPI<Scalar> {
 public:
  StepperExponential_EPI2()
  {
    this->setStepperName("EPI2");
    this->setStepperType("EPI2");
    this->setUseFSAL(false);
    this->setICConsistency("Consistent");
    this->setICConsistencyCheck(false);
    this->setZeroInitialGuess(false);
    this->setAppAction(Teuchos::null);
    this->setDefaultSolver();
    this->setOrder(2.0);
  }

  StepperExponential_EPI2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
  {
    this->setStepperName("EPI2");
    this->setStepperType("EPI2");
    this->setUseFSAL(useFSAL);
    this->setICConsistency(ICConsistency);
    this->setICConsistencyCheck(ICConsistencyCheck);
    this->setZeroInitialGuess(zeroInitialGuess);
    this->setOrder(2.0);
    this->setupTableau();

    this->setAppAction(stepperEPIAppAction);
    this->setSolver(solver);

    if (appModel != Teuchos::null) {
      this->setModel(appModel);
      this->initialize();
    }
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n";
    return Description.str();
  }

 protected:
  void setupTableau() {}
};

/// Nonmember constructor for EPI2
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperExponential_EPI2<Scalar> >
createStepperExponential_EPI2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponential_EPI2<Scalar>());
  stepper->setStepperImplicitValues(pl);
  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}


// ----------------------------------------------------------------------------
/** \brief EPI3 Definition.
 *
 *  See Tempus_StepperEPI for additional details.
 */
template <class Scalar>
class StepperExponential_EPI3 : virtual public StepperEPI<Scalar> {
 public:
  StepperExponential_EPI3()
  {
    this->setStepperName("EPI3");
    this->setStepperType("EPI3");
    this->setUseFSAL(false);
    this->setICConsistency("Consistent");
    this->setICConsistencyCheck(false);
    this->setZeroInitialGuess(false);
    this->setAppAction(Teuchos::null);
    this->setDefaultSolver();
    this->setOrder(3.0);
  }

  StepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperEPIAppAction<Scalar> >& stepperEPIAppAction)
  {
    this->setStepperName("EPI3");
    this->setStepperType("EPI3");
    this->setUseFSAL(useFSAL);
    this->setICConsistency(ICConsistency);
    this->setICConsistencyCheck(ICConsistencyCheck);
    this->setZeroInitialGuess(zeroInitialGuess);
    this->setOrder(3.0);
    this->setupTableau();

    this->setAppAction(stepperEPIAppAction);
    this->setSolver(solver);

    if (appModel != Teuchos::null) {
      this->setModel(appModel);
      this->initialize();
    }
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n";
    return Description.str();
  }

 protected:
  void setupTableau() {}
};

/// Nonmember constructor for EPI3
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperExponential_EPI3<Scalar> >
createStepperExponential_EPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperExponential_EPI3<Scalar>());
  stepper->setStepperImplicitValues(pl);
  stepper->setStepperExponentialValues(pl);

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

} // namespace Tempus

#endif // Tempus_StepperEPI_decl_hpp
