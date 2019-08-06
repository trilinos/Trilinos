// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2_decl_hpp
#define Tempus_StepperBDF2_decl_hpp

#include "Tempus_StepperImplicit.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"
#include "Tempus_StepperBDF2Observer.hpp"


namespace Tempus {

/** \brief BDF2 (Backward-Difference-Formula-2) time stepper.
 *
 *  For the implicit ODE system, \f$f(\dot{x},x,t) = 0\f$,
 *  the solution, \f$\dot{x}\f$ and \f$x\f$, is determined using a
 *  solver (e.g., a non-linear solver, like NOX).  This stepper allows
 *  for a variable time-step, \f$\Delta t\f$.  It is a 2-step method.
 *
 *  <b> Algorithm </b>
 *   - For \f$n=0\f$, set the initial condition, \f$x_0\f$.
 *   - For \f$n=1\f$, use a one-step startup stepper, e.g., Backward Euler
 *     or RK4.  The default startup stepper is 'IRK 1 Stage Theta Method'
 *     which second order.
 *   - For \f$n>1\f$, solve for \f$x_n\f$ via
 *       \f$ f\left(x_n, \dot{x}_n, t_n\right) = 0\f$
 *  where \f$
 *    \dot{x}_{n} = \frac{2\tau_n + \tau_{n-1}}{\tau_n + \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right], \f$
 *  and \f$\Delta t_n = \tau_n = t_n - t_{n-1}\f$.
 *   - \f$\dot{x}_n \leftarrow
 *    \dot{x}_{n} = \frac{2\tau_n + \tau_{n-1}}{\tau_n + \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right], \f$
 *
 *  The First-Step-As-Last (FSAL) principle is not needed BDF2.
 *  The default is to set useFSAL=false, however useFSAL=true will also work
 *  but have no affect (i.e., no-op).
 */
template<class Scalar>
class StepperBDF2 : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Constructs with a default ParameterList.
   *  - Can reset ParameterList with setParameterList().
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
  */
  StepperBDF2();

  /// Constructor
  StepperBDF2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    /// Set the stepper to use in first step
    void setStartUpStepper(std::string startupStepperName);
    void setStartUpStepper(Teuchos::RCP<Teuchos::ParameterList>startUpStepperPL=Teuchos::null);

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
    virtual Scalar getOrder() const {return order_;}
    virtual Scalar getOrderMin() const {return 1.0;}
    virtual Scalar getOrderMax() const {return 2.0;}

    virtual bool isExplicit()         const {return false;}
    virtual bool isImplicit()         const {return true;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return false;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
  //@}

  /// Return alpha = d(xDot)/dx.
  virtual Scalar getAlpha(const Scalar dt) const {return getAlpha(dt,dt);}
  virtual Scalar getAlpha(const Scalar dt, const Scalar dtOld) const
    { return (Scalar(2.0)*dt + dtOld)/(dt*(dt + dtOld)); }
  /// Return beta  = d(x)/dx.
  virtual Scalar getBeta (const Scalar   ) const { return Scalar(1.0); }

  /// Compute the first time step given the supplied startup stepper
  virtual void computeStartUp(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

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

  Teuchos::RCP<Stepper<Scalar> >             startUpStepper_;
  Teuchos::RCP<StepperBDF2Observer<Scalar> > stepperBDF2Observer_;
  Scalar                                     order_;
};

/** \brief Time-derivative interface for BDF2.
 *
 *  Given the state \f$x_n\f$, compute the BDF2 time-derivative,
 *  \f[
 *    \dot{x}_{n} = \frac{2\tau_n + \tau_{n-1}}{\tau_n + \tau_{n-1}}
 *                  \left[ \frac{x_n-x_{n-1}}{\tau_n}\right]
 *                -  \frac{\tau_n}{\tau_n + \tau_{n-1}}
 *                   \left[ \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right]
 *  \f]
 *  where
 *  \f[
 *   \tau_n = t_n - t_{n-1}.
 *   \f]
 *  \f$\ddot{x}\f$ is not used and set to null.
 */
template <typename Scalar>
class StepperBDF2TimeDerivative
  : virtual public Tempus::TimeDerivative<Scalar>
{
public:

  /// Constructor
  StepperBDF2TimeDerivative(
    Scalar dt, Scalar dtOld, Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOldOld)
  { initialize(dt, dtOld, xOld, xOldOld); }

  /// Destructor
  virtual ~StepperBDF2TimeDerivative() {}

  /// Compute the time derivative.
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null)
  {
    xDotDot = Teuchos::null;
    // Calculate the BDF2 x dot vector
    const Scalar a = ((Scalar(2.0)*dt_ + dtOld_)/(dt_ + dtOld_))/dt_;
    const Scalar b = (                       dt_/(dt_ + dtOld_))/dtOld_;
    //xDot = a*(x_n - x_{n-1}) - b*(x_{n-1} - x_{n-2})
    Thyra::V_StVpStV(xDot.ptr(), a, *x, -(a+b), *xOld_);
    Thyra::Vp_StV(xDot.ptr(), b, *xOldOld_);
  }

  virtual void initialize(Scalar dt, Scalar dtOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld,
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOldOld)
  { dt_ = dt; dtOld_ = dtOld; xOld_ = xOld; xOldOld_ = xOldOld;}

private:

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOld_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xOldOld_;
  Scalar                                         dt_;    // = t_n - t_{n-1}
  Scalar                                         dtOld_; // = t_{n-1} - t_{n-2}
};


} // namespace Tempus

#endif // Tempus_StepperBDF2_decl_hpp
