// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBackwardEuler_decl_hpp
#define Tempus_StepperBackwardEuler_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperBackwardEulerObserver.hpp"
#include "Tempus_StepperOptimizationInterface.hpp"


namespace Tempus {

/** \brief Backward Euler time stepper.
 *
 *  For the implicit ODE system, \f$\mathcal{F}(\dot{x},x,t) = 0\f$,
 *  the solution, \f$\dot{x}\f$ and \f$x\f$, is determined using a
 *  solver (e.g., a non-linear solver, like NOX).
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Backward Euler is simply,
 *   - Solve \f$f(\dot{x}=(x_n-x_{n-1})/\Delta t_n, x_n, t_n)=0\f$ for \f$x_n\f$
 *   - \f$\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n\f$
 *
 *  The First-Step-As-Last (FSAL) principle is not needed with Backward Euler.
 *  The default is to set useFSAL=false, however useFSAL=true will also work
 *  but have no affect (i.e., no-op).
 */
template<class Scalar>
class StepperBackwardEuler :
    virtual public Tempus::StepperImplicit<Scalar>,
    virtual public Tempus::StepperOptimizationInterface<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
  */
  StepperBackwardEuler();

  /// Constructor
  StepperBackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperObserver<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual Teuchos::RCP<StepperObserver<Scalar> > getObserver() const
    { return this->stepperBEObserver_; }

    /// Set the predictor
    void setPredictor(std::string predictorType = "None");
    void setPredictor(Teuchos::RCP<Stepper<Scalar> > predictorStepper);

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
    virtual Scalar getOrder() const {return 1.0;}
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 1.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
  //@}

    /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const { return Scalar(1.0)/dt; }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar   ) const { return Scalar(1.0); }

  /// Compute predictor given the supplied stepper
  virtual void computePredictor(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Implementation of StepperOptimizationInterface
  //@{
    virtual int stencilLength() const;
    virtual void computeStepResidual(
      Thyra::VectorBase<Scalar>& residual,
      const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t,
      const Thyra::VectorBase<Scalar>& p,
      const int param_index) const;
    virtual void computeStepJacobian(
      Thyra::LinearOpBase<Scalar>& jacobian,
      const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t,
      const Thyra::VectorBase<Scalar>& p,
      const int param_index,
      const int deriv_index) const;
    virtual void computeStepParamDeriv(
      Thyra::LinearOpBase<Scalar>& deriv,
      const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t,
      const Thyra::VectorBase<Scalar>& p,
      const int param_index) const;
    virtual void computeStepSolver(
      Thyra::LinearOpWithSolveBase<Scalar>& jacobian_solver,
      const Teuchos::Array< Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x,
      const Teuchos::Array<Scalar>& t,
      const Thyra::VectorBase<Scalar>& p,
      const int param_index) const;
  //@}

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

  Teuchos::RCP<Stepper<Scalar> >                      predictorStepper_;
  Teuchos::RCP<StepperBackwardEulerObserver<Scalar> > stepperBEObserver_;

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
class StepperBackwardEulerTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperBackwardEulerTimeDerivative(
    Scalar s, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld)
  { initialize(s, xOld); }

  /// Destructor
  virtual ~StepperBackwardEulerTimeDerivative() {}

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


} // namespace Tempus

#endif // Tempus_StepperBackwardEuler_decl_hpp
