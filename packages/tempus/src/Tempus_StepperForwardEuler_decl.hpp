// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperForwardEuler_decl_hpp
#define Tempus_StepperForwardEuler_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperExplicit.hpp"
#include "Tempus_StepperForwardEulerObserver.hpp"


namespace Tempus {

/** \brief Forward Euler time stepper.
 *
 *  For the explicit ODE system,
 *  \f[
 *    \dot{x} = \bar{f}(x,t),
 *  \f]
 *  the Forward Euler stepper can be written as
 *  \f[
 *    x_{n} = x_{n-1} + \Delta t\, \bar{f}(x_{n-1},t_{n-1})
 *  \f]
 *  Forward Euler is an explicit time stepper (i.e., no solver used).
 *  Note that the time derivative by definition is
 *  \f[
 *    \dot{x}_{n} = \bar{f}(x_{n},t_{n}),
 *  \f]
 *
 *  <b> Algorithm </b>
 *  The single-timestep algorithm for Forward Euler is simply,
 *   - \f$\dot{x}_{n-1} \leftarrow \bar{f}(x_{n-1},t_{n-1})\f$
 *   - \f$x_{n} \leftarrow x_{n-1} + \Delta t\, \dot{x}_{n-1}\f$
 *
 *  Note that \f$x_n\f$ and \f$\dot{x}_{n-1}\f$ are not at the same time
 *  level at the end of the time step (i.e., they are not sync'ed).
 *
 *  To have them at the same time level, we can use the First-Step-As-Last
 *  (FSAL) principle where the function evaulation from the last time step
 *  can be used as the first function evalulation of the current step.
 *  For the Forward Euler, the FSAL algorithm is
 *   - \f$x_{n} \leftarrow x_{n-1} + \Delta t\, \dot{x}_{n-1}\f$
 *   - \f$\dot{x}_n \leftarrow \bar{f}(x_{n},t_{n})\f$
 *
 *  The default for Forward Euler is to use FSAL (useFSAL=true).
 */
template<class Scalar>
class StepperForwardEuler : virtual public Tempus::StepperExplicit<Scalar>
{
public:

  /** \brief Default constructor.
   *
   *  - Constructs with a default ParameterList.
   *  - Can reset ParameterList with setParameterList().
   *  - Requires subsequent setModel() and initialize() calls before calling
   *    takeStep().
  */
  StepperForwardEuler();

  /// Constructor
  StepperForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null);

  /// \name Basic stepper methods
  //@{
    virtual void setObserver(
      Teuchos::RCP<StepperObserver<Scalar> > obs = Teuchos::null);

    /// Initialize during construction and after changing input parameters.
    virtual void initialize();

    /// Set the initial conditions, make them consistent, and set needed memory.
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

    virtual OrderODE getOrderODE()   const {return FIRST_ORDER_ODE;}
  //@}

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

protected:

  Teuchos::RCP<StepperForwardEulerObserver<Scalar> >  stepperFEObserver_;

};

} // namespace Tempus

#endif // Tempus_StepperForwardEuler_decl_hpp
