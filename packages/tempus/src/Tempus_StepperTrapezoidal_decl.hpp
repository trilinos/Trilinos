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
#include "Tempus_StepperTrapezoidalObserver.hpp"


namespace Tempus {

/** \brief Trapezoidal method time stepper.
 *
 *  For the implicit ODE system, \f$\mathcal{F}(\dot{x},x,t) = 0\f$,
 *  the solution, \f$\dot{x}\f$ and \f$x\f$, is determined using a
 *  solver (e.g., a non-linear solver, like NOX).
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Trapezoidal method is simply,
 *   - Solve \f$f(\dot{x}=(x_n-x_{n-1})/(\Delta t_n/2) - \dot{x}_{n-1}, x_n, t_n)=0\f$ for \f$x_n\f$
 *   - \f$\dot{x}_n \leftarrow (x_n-x_{n-1})/(\Delta t_n/2) - \dot{x}_{n-1}\f$
 *
 *   The First-Step-As-Last (FSAL) principle is required for the Trapezoidal
 *   Stepper (i.e., useFSAL=true)!  There are at least two ways around this,
 *   but are not implemented.
 *    - Do a solve for xDotOld, xDot_{n-1}, at each time step as for the
 *      initial conditions.  This is expensive since you would be doing
 *      two solves every time step.
 *    - Use evaluateExplicitODE to get xDot_{n-1} if the application
 *      provides it.  Explicit evaluations are cheaper but requires the
 *      application to implement xDot = f(x,t).
 */
template<class Scalar>
class StepperTrapezoidal : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Constructs with a default ParameterList.
   *  - Can reset ParameterList with setParameterList().
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
  */
  StepperTrapezoidal();

  /// Constructor
  StepperTrapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

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
    virtual Scalar getOrder() const {return 2.0;}
    virtual Scalar getOrderMin() const {return 2.0;}
    virtual Scalar getOrderMax() const {return 2.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
  //@}

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const { return Scalar(2.0)/dt; }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar   ) const { return Scalar(1.0); }

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

private:

  Teuchos::RCP<Stepper<Scalar> >                    predictorStepper_;
  Teuchos::RCP<StepperTrapezoidalObserver<Scalar> > stepperTrapObserver_;

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
