// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperLeapfrog_decl_hpp
#define Tempus_StepperLeapfrog_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_StepperLeapfrogObserver.hpp"


namespace Tempus {

/** \brief Leapfrog time stepper.
 *
 *  For the governing equation,
 *  \f[
 *    M(t) \ddot{x} + K(t) x = F(t),
 *  \f]
 *  one can write the explicit ODE system,
 *  \f[
 *    \ddot{x} = f(x,t),
 *  \f]
 *  where
 *  \f[
 *    f(x,t) = \left(M(t)\right)^{-1} \left( F(t) - K(t) x \right).
 *  \f]
 *  The Leapfrog stepper can be written as
 *  \f{eqnarray*}{
 *    x_{n+1}         & = & x_{n} + \Delta t\, \dot{x}_{n+1/2} \\
 *    \ddot{x}_{n+1}  & = & f(x_{n+1},t_{n+1}) \\
 *    \dot{x}_{n+3/2} & = & \dot{x}_{n+1/2} + \Delta t\, \ddot{x}_{n+1}
 *  \f}
 *  where the position and velocity are leapfrogged over each other.
 *  On startup the velocity half-step can be obtained with
 *  \f{eqnarray*}{
 *    \dot{x}_{n+1/2} & = & \dot{x}_{n} + \frac{1}{2} \Delta t\, \ddot{x}_{n} \\
 *    \dot{x}_{n+1/2} & = & \dot{x}_{n} + \frac{1}{2} \Delta t\, f(x_{n},t_{n})
 *  \f}
 *  and to complete the time step, the final velocity half-step is obtained
 *  with
 *  \f{eqnarray*}{
 *    \dot{x}_{n+1}&=&\dot{x}_{n+1/2} +\frac{1}{2} \Delta t\, \ddot{x}_{n+1} \\
 *    \dot{x}_{n+1}&=&\dot{x}_{n+1/2} +\frac{1}{2} \Delta t\, f(x_{n+1},t_{n+1})
 *  \f}
 *
 *  <b> Algorithm </b>
 *
 *  Beginning with \f$(x_n,\dot{x}_{n+1/2})\f$ or \f$(x_n,\dot{x}_{n})\f$
 *  and/or ending with \f$(x_{n+1},\dot{x}_{n+3/2})\f$ or
 *  \f$(x_{n+1},\dot{x}_{n+1})\f$, the algorithm for Leapfrog is
 *   - if "startup"
 *     - \f$ \ddot{x}_{n} \leftarrow f(x_{n},t_{n}) \f$
 *     - \f$ \dot{x}_{n+1/2} \leftarrow
 *           \dot{x}_{n} + \frac{1}{2} \Delta t\, \ddot{x}_{n} \f$
 *   - \f$x_{n+1} \leftarrow x_{n} + \Delta t\, \dot{x}_{n+1/2} \f$
 *   - \f$\ddot{x}_{n+1} \leftarrow f(x_{n+1},t_{n+1}) \f$
 *   - if "ending"
 *     - \f$ \dot{x}_{n+1} \leftarrow
 *           \dot{x}_{n+1/2} +\frac{1}{2} \Delta t\, \ddot{x}_{n+1} \f$
 *   - else
 *     - \f$ \dot{x}_{n+3/2} \leftarrow
 *           \dot{x}_{n+1/2} + \Delta t\, \ddot{x}_{n+1} \f$
 */
template<class Scalar>
class StepperLeapfrog : virtual public Tempus::Stepper<Scalar>
{
public:

  /// Constructor
  StepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel);
    virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
      getModel(){return appModel_;}

    virtual void setSolver(std::string solverName);
    virtual void setSolver(
      Teuchos::RCP<Teuchos::ParameterList> solverPL=Teuchos::null);
    virtual void setSolver(
        Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver);
    virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > getSolver() const
      { return Teuchos::null; }
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    /// Pass initial guess to Newton solver (only relevant for implicit solvers)
    virtual void setInitialGuess(Teuchos::RCP<const Thyra::VectorBase<Scalar> > initial_guess)
       {initial_guess_ = initial_guess;}

    virtual std::string getStepperType() const
     { return stepperPL_->get<std::string>("Stepper Type"); }

    /// Get a default (initial) StepperState
    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
    virtual Scalar getOrder() const {return 2.0;}
    virtual Scalar getOrderMin() const {return 2.0;}
    virtual Scalar getOrderMax() const {return 2.0;}
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) const
      {return Scalar(1.0e+99);}

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return false;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}
  //@}

  virtual void setIsXDotXDotInitialized(bool tf)
  { stepperPL_->set<bool>("Is xDotDot Initialized", int(tf)); }
  virtual bool getIsXDotXDotInitialized() const
  { return stepperPL_->get<bool>("Is xDotDot Initialized"); }

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
  StepperLeapfrog();

protected:

  Teuchos::RCP<Teuchos::ParameterList>               stepperPL_;
  /// Explicit ODE ModelEvaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appModel_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>          inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>         outArgs_;

  Teuchos::RCP<StepperObserver<Scalar> >             stepperObserver_;
  Teuchos::RCP<StepperLeapfrogObserver<Scalar> >     stepperLFObserver_;

  Teuchos::RCP<const Thyra::VectorBase<Scalar> >      initial_guess_;

};

} // namespace Tempus

#endif // Tempus_StepperLeapfrog_decl_hpp
