#ifndef TEMPUS_STEPPER_HPP
#define TEMPUS_STEPPER_HPP

//Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
// Thyra
#include "Thyra_ModelEvaluator.hpp"
// Tempus
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {


/** \brief Thyra Base interface for time steppers.
 *
 * <b>Design Considerations</b>
 *   - Time steppers are designed to take a single "basic" time step, and
 *     thus everything refers to the current time step.
 *   - A "basic" time step can refer to
 *     - a single explicit time step
 *     - a single implicit solve for a time step
 *     - a single solve for a IMEX time step
 *   - Multiple time steps should be managed by Integrators.
 *   - Steppers can be built from other Steppers.
 *     - An operator-split Stepper is possible with interoperable Steppers.
 *   - For explicit steppers, only one ModelEvaluator and one solution
 *     vector are required.
 *   - For implicit steppers, only one ModelEvaluator, one solution
 *     vector, and one solver are required.
 *   - For more complex steppers, multiple ModelEvaluators, solution
 *     vectors, and solvers are possible when a common single time-integration
 *     method is desired for all solutions. Examples:
 *     - Solution A with ModelEvaluator A and Solution B with ModelEvaluator B
 *       using the same solver
 *     - Solution A with ModelEvaluator A using Solver A and Solution B with
 *       ModelEvaluator B using Solver B
 *     - Solution A with ModelEvaluator A using Solver A and Solutions A and B
 *       with ModelEvaluator C using Solver B
 *   - Steppers may maintain their own time history of the solution, e.g.,
 *     BDF steppers.
 */
template<class Scalar>
class Stepper
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<Stepper<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /// Destructor
//  virtual ~Stepper();

  /// \name Basic stepper methods
  //@{
    /// Take the specified timestep, dt, and return true if successful.
    virtual bool takeStep(const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    virtual void setStepperState(
      const Teuchos::RCP<Tempus::StepperState<Scalar> >& stepperState) = 0;

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getStepperState() = 0;
  //@}

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
//    virtual void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pl);
//    virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
//    virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
//    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
//    virtual std::string description() const;
//    virtual void describe(Teuchos::FancyOStream        & out,
//                          const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Error estimation methods
  //@{

  //@}

  /// \name Observer methods
  //@{

  //@}

  /// \name Adjoint methods
  //@{
  //virtual Scalar takeAdjointStep();
  //@}

  /// \name Solution history methods
  //@{

  /// Functionality like InterpolationBuffer for multi-step methods, BDF.

  //@}

};
} // namespace Tempus
#endif // TEMPUS_STEPPER_HPP
