#ifndef TEMPUS_STEPPERBACKWARDEULER_HPP
#define TEMPUS_STEPPERBACKWARDEULER_HPP

#include "Tempus_Stepper.hpp"
#include "Tempus_ResidualModelEvaluator.hpp"
#include "NOX_Thyra.H"

namespace Tempus {


/** \brief Backward Euler time stepper.
 *  Backward Euler is an implicit time stepper (i.e., with solver).
 */
template<class Scalar>
class StepperBackwardEuler : virtual public Tempus::Stepper<Scalar>
{
public:

  /// Constructor
  StepperBackwardEuler(
    Teuchos::RCP<Teuchos::ParameterList>                pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel );

  /// \name Basic stepper methods
  //@{
    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory);

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState();
  //@}

  /// \name ParameterList methods
  //@{
    void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
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

  Teuchos::RCP<Teuchos::ParameterList>              pList_;
  Teuchos::RCP<ResidualModelEvaluator<Scalar> >     residualModel_;
  //Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  Teuchos::RCP<Thyra::NOXNonlinearSolver>  solver_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar>  inArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> outArgs_;

  // Compute the balancing time derivative as a function of x
  std::function<void (const Thyra::VectorBase<Scalar> &,
                            Thyra::VectorBase<Scalar> &)>
  xDotFunction(Scalar dt,Teuchos::RCP<const Thyra::VectorBase<Scalar> > x_old);

};
} // namespace Tempus

#include "Tempus_StepperBackwardEuler_impl.hpp"

#endif // TEMPUS_STEPPERBACKWARDEULER_HPP
