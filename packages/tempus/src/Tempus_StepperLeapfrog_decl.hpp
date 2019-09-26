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
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperObserverComposite.hpp"
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
 *
 *  The First-Step-As-Last (FSAL) principle is not used with Leapfrog
 *  because of the algorithm's prescribed order of solution update.
 *  The default is to set useFSAL=false, however useFSAL=true will also
 *  work (i.e., no-op), but issue a warning that it will have no affect.
 */
template<class Scalar>
class StepperLeapfrog : virtual public Tempus::StepperExplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
  */
  StepperLeapfrog();

  /// Constructor
  StepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperObserver<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    virtual Teuchos::RCP<StepperObserver<Scalar> > getObserver() const
    { return this->stepperObserver_; }

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
    virtual Scalar getInitTimeStep(
        const Teuchos::RCP<SolutionHistory<Scalar> >& /* solutionHistory */) const
      {return Scalar(1.0e+99);}

    virtual bool isExplicit()         const {return true;}
    virtual bool isImplicit()         const {return false;}
    virtual bool isExplicitImplicit() const
      {return isExplicit() and isImplicit();}
    virtual bool isOneStepMethod()   const {return true;}
    virtual bool isMultiStepMethod() const {return !isOneStepMethod();}

    virtual OrderODE getOrderODE()   const {return SECOND_ORDER_ODE;}
  //@}

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  std::string getICConsistencyDefault() const { return "Consistent"; }

  /// \name Overridden from Teuchos::Describable
  //@{
    virtual void describe(Teuchos::FancyOStream        & out,
                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

protected:

  Teuchos::RCP<StepperObserverComposite<Scalar> >    stepperObserver_;
  Teuchos::RCP<StepperLeapfrogObserver<Scalar> >     stepperLFObserver_;

};

} // namespace Tempus

#endif // Tempus_StepperLeapfrog_decl_hpp
