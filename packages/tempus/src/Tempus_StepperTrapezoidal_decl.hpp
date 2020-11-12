// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperTrapezoidal_decl_hpp
#define Tempus_StepperTrapezoidal_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperTrapezoidalAppAction.hpp"
#include "Tempus_StepperOptimizationInterface.hpp"


namespace Tempus {

/** \brief Trapezoidal method time stepper.
 *
 *  For the implicit ODE system, \f$\mathcal{F}(\dot{x},x,t) = 0\f$,
 *  the solution, \f$\dot{x}\f$ and \f$x\f$, is determined using a
 *  solver (e.g., a non-linear solver, like NOX).
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Forward Euler is
 *  \f{algorithm}{
 *  \renewcommand{\thealgorithm}{}
 *  \caption{Forward Euler}
 *  \begin{algorithmic}[1]
 *    \State {\it appAction.execute(solutionHistory, stepper, BEGIN\_STEP)}
 *    \State {Get $\dot{x}$ from SolutionHistory or from Stepper}
 *    \State {\it appAction.execute(solutionHistory, stepper, BEFORE\_SOLVE)}
 *    \State $\mathcal{F}_n(\dot{x}=(x_n-x_{n-1})/(\Delta t_n/2) - \dot{x}_{n-1}, x_n, t_n)=0\f$ for \f$x_n$
 *    \State {\it appAction.execute(solutionHistory, stepper, AFTER\_SOLVE)}
 *    \State $\dot{x}_n \leftarrow (x_n-x_{n-1})/(\Delta t_n/2) - \dot{x}_{n-1}$
 *    \State {\it appAction.execute(solutionHistory, stepper, END\_STEP)}
 *  \end{algorithmic}
 *  \f}
 *
 *  The First-Same-As-Last (FSAL) principle is required for the Trapezoidal
 *  Stepper (i.e., useFSAL=true)!  With useFSAL=true does assume that the
 *  solution, \f$x\f$, and its time derivative, \f$\dot{x}\f$, are consistent
 *  at the initial conditions (ICs), i.e.,
 *  \f$\dot{x}_{0} = \bar{f}(x_{0},t_{0})\f$.  This can be ensured by setting
 *  setICConsistency("Consistent"), and checked with
 *  setICConsistencyCheck(true).
 *
 *  <b> Iteration Matrix, \f$W\f$.</b>
 *  Recalling that the definition of the iteration matrix, \f$W\f$, is
 *  \f[
 *    W = \alpha \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \beta  \frac{\partial \mathcal{F}_n}{\partial x_n},
 *  \f]
 *  where \f$ \alpha \equiv \frac{\partial \dot{x}_n(x_n) }{\partial x_n}, \f$
 *  and \f$ \beta \equiv \frac{\partial x_n}{\partial x_n} = 1\f$, and
 *  the time derivative for Trapezoidal method is
 *  \f[
 *    \dot{x}_n = (x_n-x_{n-1})/(\Delta t/2) - \dot{x}_{n-1},
 *  \f]
 *  we can determine that
 *  \f$ \alpha = \frac{2}{\Delta t} \f$
 *  and \f$ \beta = 1 \f$, and therefore write
 *  \f[
 *    W = \frac{2}{\Delta t}
 *        \frac{\partial \mathcal{F}_n}{\partial \dot{x}_n}
 *      + \frac{\partial \mathcal{F}_n}{\partial x_n}.
 *  \f]
 */
template<class Scalar>
class StepperTrapezoidal : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  Requires subsequent setModel(), setSolver() and initialize()
   *  calls before calling takeStep().
  */
  StepperTrapezoidal();

  /// Constructor
  StepperTrapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> >& stepperTrapAppAction);

  /// \name Basic stepper methods
  //@{
    virtual void setAppAction(
    Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> > appAction);

    virtual Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> > getAppAction() const
    { return stepperTrapAppAction_; }

    /// Set the initial conditions and make them consistent.
    virtual void setInitialConditions (
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {return 2.0;}
    virtual Scalar getOrderMin() const {return 2.0;}
    virtual Scalar getOrderMax() const {return 2.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
    virtual void setUseFSAL(bool a) { this->setUseFSALTrueOnly(a); }
    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
  //@}

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const { return Scalar(2.0)/dt; }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar   ) const { return Scalar(1.0); }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  virtual bool isValidSetup(Teuchos::FancyOStream & out) const;


private:

  Teuchos::RCP<StepperTrapezoidalAppAction<Scalar> > stepperTrapAppAction_;

};

/** \brief Time-derivative interface for Trapezoidal method.
 *
 *  Given the state \f$x\f$, compute the Trapezoidal method time-derivative,
 *  \f[
 *    \dot{x}_{n} = \frac{(x_{n} - x_{n-1})}{(\Delta t_n/2)} - \dot{x}_{n-1}.
 *  \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperTrapezoidalTimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperTrapezoidalTimeDerivative( Scalar s,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotOld)
  { initialize(s, xOld, xDotOld); }

  /// Destructor
  virtual ~StepperTrapezoidalTimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;
    // Calculate the Trapezoidal method x dot vector
    Thyra::V_StVpStV(xDot.ptr(),s_,*x,-s_,*xOld_);
    Thyra::V_VpStV  (xDot.ptr(),*xDot,Scalar(-1.0),*xDotOld_);
  }

  virtual void initialize(Scalar s,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotOld)
  { s_ = s; xOld_ = xOld; xDotOld_ = xDotOld; }

private:

  Scalar                                         s_;    // = 1.0/(dt/2)
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xDotOld_;
};


} // namespace Tempus

#endif // Tempus_StepperTrapezoidal_decl_hpp
