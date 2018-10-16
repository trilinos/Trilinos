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
 *  for a variable time-step \f$dt\f$.  It is a 3-step method.
 *
 *  <b> Algorithm </b>
 *   - Select initial guess \f$x_n\f$ for \f$n\f$.
 *   - Compute \f$x_n\f$ for n=1 using some time-integration scheme,
 *     e.g., Backward Euler or RK4.
 *   - Solve
 *   \f[
 *   f\left(\dot{x} = \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}
 *                 + \left(\frac{1}{\tau_{n-1}-\tau_n} \right)
 *                   \left(\frac{x_n-x_{n-1}}{\tau_{n-1} + \tau_n}
 *                         - \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right)
 *                  (2\tau_n + \tau_{n-1}), x_n, t_n\right) = 0
 *   \f]
 *   for \f$x_n\f$ (n > 1) where
 *   \f[
 *   \tau_n = t_n - t_{n-1}.
 *   \f]
 */
template<class Scalar>
class StepperBDF2 : virtual public Tempus::StepperImplicit<Scalar>
{
public:

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
  //@}

  /// Compute the first time step given the supplied startup stepper
  virtual void computeStartUp(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);
    
  /// Pass initial guess to Newton solver 
  virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
       {initial_guess_ = initial_guess;}


  /// Provide temporary xDot memory for Stepper if SolutionState doesn't.
  virtual Teuchos::RCP<Thyra::VectorBase<Scalar> > getXDotTemp(
    Teuchos::RCP<Thyra::VectorBase<Scalar> > x);

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

  /// Default Constructor -- not allowed
  StepperBDF2();

private:

  Teuchos::RCP<Stepper<Scalar> >             startUpStepper_;

  Teuchos::RCP<StepperObserver<Scalar> >     stepperObserver_;
  Teuchos::RCP<StepperBDF2Observer<Scalar> > stepperBDF2Observer_;
  Scalar                                     order_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >   xDotTemp_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;  
};

/** \brief Time-derivative interface for BDF2.
 *
 *  Given the state \f$x\f$, compute the BDF2 time-derivative,
 *  \f[
 *    \dot{x}_{n} = \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}
 *                + \left(\frac{1}{\tau_{n-1}-\tau_n} \right)
 *                  \left(\frac{x_n-x_{n-1}}{\tau_{n-1} + \tau_n}
 *                        - \frac{x_{n-1}-x_{n-2}}{\tau_{n-1}}\right)
 *                  (2\tau_n + \tau_{n-1})
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
    const Scalar a = (1.0/(dt_ + dtOld_))*(2.0*dt_ + dtOld_)/dt_;
    const Scalar b = (1.0/(dt_ + dtOld_))*(dt_/dtOld_);
    //xDot = a*(x_n-x_{n-1})
    Thyra::V_StVpStV(xDot.ptr(),a,*x,-a,*xOld_);
    Teuchos::RCP<Thyra::VectorBase<Scalar> > tmp =
      Thyra::createMember<Scalar>(x->space());
    //tmp = b*(x_{n-1} - x_{n-2})
    Thyra::V_StVpStV(tmp.ptr(),b,*xOld_,-b,*xOldOld_);
    //xDot = xDot - tmp;
    Thyra::Vp_StV(xDot.ptr(), -1.0, *tmp);
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
