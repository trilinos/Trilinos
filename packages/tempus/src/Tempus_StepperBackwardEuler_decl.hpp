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
 *   - \f$\dot{x}_n \leftarrow (x_n-x_{n-1})/\Delta t_n\f$ [Optional]
 *   - Solve \f$f(\dot{x}_n,x_n,t_n)=0\f$ for \f$\dot{x}_n\f$ [Optional]
 */
template<class Scalar>
class StepperBackwardEuler : virtual public Tempus::StepperImplicit<Scalar>
{
public:

  /// Constructor
  StepperBackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    /// Set the predictor
    void setPredictor(std::string predictorName);
    void setPredictor(Teuchos::RCP<Teuchos::ParameterList>predPL=Teuchos::null);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

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
  //@}

  /// Compute predictor given the supplied stepper
  virtual void computePredictor(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

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
  StepperBackwardEuler();

private:

  Teuchos::RCP<Stepper<Scalar> >                      predictorStepper_;
  Teuchos::RCP<StepperObserver<Scalar> >              stepperObserver_;
  Teuchos::RCP<StepperBackwardEulerObserver<Scalar> > stepperBEObserver_;

  Teuchos::RCP<Thyra::VectorBase<Scalar> >            xDotTemp_;
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
