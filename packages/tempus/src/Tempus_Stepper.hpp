#ifndef Tempus_Stepper_hpp
#define Tempus_Stepper_hpp

//Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
// Thyra
#include "Thyra_ModelEvaluator.hpp"
// Tempus
#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"


namespace Tempus {


/** \brief Thyra Base interface for time steppers.
 *
 * <b>Design Considerations</b>
 *   - Time steppers are designed to take a single time step.
 *     - a single implicit solve for a time step
 *     - a single solve for a IMEX time step
 *   - Multiple time steps should be managed by Integrators.
 *   - Steppers can be built from other Sub-Steppers.
 *     - An operator-split Stepper is possible with interoperable Steppers.
 *   - For explicit steppers, only one ModelEvaluator and one solution
 *     vector are required.
 *   - For implicit steppers, only one ModelEvaluator, one solution
 *     vector, and one solver are required.
 *   - Steppers will PASS/FAIL the time step based on Solver, error and
 *     order requirements, and not adjust the time step size.
 *   - Steppers can provide a suggested time step size for the next time step.
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
 *
 * <b> CS Design Considerations</b>
 *   - Steppers will be fully constructed with input or default parameters.
 *   - All input parameters (i.e., ParameterList) can be set by public methods.
 *   - The Stepper ParameterList must be consistent.
 *     - The "set" methods which update parameters in the ParameterList
 *       must update the Stepper ParameterList.
 */
template<class Scalar>
class Stepper
  : virtual public Teuchos::Describable,
    virtual public Teuchos::VerboseObject<Stepper<Scalar> >,
    virtual public Teuchos::ParameterListAcceptor
{
public:

  /// \name Basic stepper methods
  //@{
    virtual void setModel(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel
      ) = 0;
    virtual void setNonConstModel(
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& transientModel) = 0;

    /// Initialize during construction and after changing input parameters.
    virtual void initialize() = 0;

    /// Take the specified timestep, dt, and return true if successful.
    virtual void takeStep(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory) = 0;

    virtual Teuchos::RCP<Tempus::StepperState<Scalar> > getDefaultStepperState()
       = 0;
    virtual Scalar getOrder() const = 0;
    virtual Scalar getOrderMin() const = 0;
    virtual Scalar getOrderMax() const = 0;
  //@}

  /// \name Helper functions
  //@{
    /// Validate that the model supports explicit ODE evaluation, f(x,t) [=xdot]
    /** Currently the convention to evaluate f(x,t) is to set xdot=null!
     *  There is no InArgs support to test if xdot is null, so we set
     *  xdot=null and hopefully the ModelEvaluator can handle it.
     */
    void validExplicitODE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const
    {
      TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
      typedef Thyra::ModelEvaluatorBase MEB;
      const MEB::InArgs<Scalar>  inArgs  = model->createInArgs();
      const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
      const bool supports = inArgs.supports(MEB::IN_ARG_x) and
                            outArgs.supports(MEB::OUT_ARG_f);

      TEUCHOS_TEST_FOR_EXCEPTION( supports == false, std::logic_error,
        model->description() << "can not support an explicit ODE with\n"
        << "  IN_ARG_x  = " << inArgs.supports(MEB::IN_ARG_x) << "\n"
        << "  OUT_ARG_f = " << outArgs.supports(MEB::OUT_ARG_f) << "\n"
        << "Explicit ODE requires:\n"
        << "  IN_ARG_x  = true\n"
        << "  OUT_ARG_f = true\n"
        << "\n"
        << "NOTE: Currently the convention to evaluate f(x,t) is to set\n"
        << "xdot=null!  There is no InArgs support to test if xdot is null,\n"
        << "so we set xdot=null and hope the ModelEvaluator can handle it.\n");

      return;
    }

    /// Validate that the model supports implicit ODE/DAE evaluation, f(xdot,x,t) [= 0]
    void validImplicitODE_DAE(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model) const
    {
      TEUCHOS_TEST_FOR_EXCEPT( is_null(model) );
      typedef Thyra::ModelEvaluatorBase MEB;
      const MEB::InArgs<Scalar>  inArgs  = model->createInArgs();
      const MEB::OutArgs<Scalar> outArgs = model->createOutArgs();
      const bool supports = inArgs.supports(MEB::IN_ARG_x) and
                            inArgs.supports(MEB::IN_ARG_x_dot) and
                            inArgs.supports(MEB::IN_ARG_alpha) and
                            inArgs.supports(MEB::IN_ARG_beta) and
                            outArgs.supports(MEB::OUT_ARG_f) and
                            outArgs.supports(MEB::OUT_ARG_W);

      TEUCHOS_TEST_FOR_EXCEPTION( supports == false, std::logic_error,
        model->description() << "can not support an implicit ODE with\n"
        << "  IN_ARG_x     = " << inArgs.supports(MEB::IN_ARG_x) << "\n"
        << "  IN_ARG_x_dot = " << inArgs.supports(MEB::IN_ARG_x_dot) << "\n"
        << "  IN_ARG_alpha = " << inArgs.supports(MEB::IN_ARG_alpha) << "\n"
        << "  IN_ARG_beta  = " << inArgs.supports(MEB::IN_ARG_beta) << "\n"
        << "  OUT_ARG_f    = " << outArgs.supports(MEB::OUT_ARG_f) << "\n"
        << "  OUT_ARG_W    = " << outArgs.supports(MEB::OUT_ARG_W) << "\n"
        << "Implicit ODE requires:\n"
        << "  IN_ARG_x     = true\n"
        << "  IN_ARG_x_dot = true\n"
        << "  IN_ARG_alpha = true\n"
        << "  IN_ARG_beta  = true\n"
        << "  OUT_ARG_f    = true\n"
        << "  OUT_ARG_W    = true\n");

      return;
    }
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
#endif // Tempus_Stepper_hpp
