#ifndef TEMPUS_TYPES_HPP
#define TEMPUS_TYPES_HPP

namespace Tempus {


/** \brief Step types for integrators and steppers.
 *  Integrators primarily have CONSTANT_STEP_SIZE and VARIABLE_STEP_SIZE
 *  step types.  However steppers can have
 *  - CONSTANT_STEP_SIZE and VARIABLE_STEP_SIZE (e.g., BDF methods)
 *  - UNMODIFIABLE_STEP_SIZE, do not let the stepper change the time step given from the integrator
 *  - MODIFIABLE_STEP_SIZE, let the stepper adjust the time step to meet stepper order and error tolerances
 *  - MODIFIABLE_UP_STEP_SIZE, let the stepper increase the time step (e.g., aggressively increase time step)
 *  - MODIFIABLE_DOWN_STEP_SIZE, let the stepper decrease the time step (e.g., ensure not to miss an output time)
 */
enum StepType {
  CONSTANT_STEP_SIZE,       ///< Step size is constant.
  VARIABLE_STEP_SIZE,       ///< Step size is variable.
  UNMODIFIABLE_STEP_SIZE,   ///< Step size is not modifiable.
  MODIFIABLE_STEP_SIZE,     ///< Step size is modifiable.
  MODIFIABLE_UP_STEP_SIZE,  ///< Step size is modifiable up (increase).
  MODIFIABLE_DOWN_STEP_SIZE ///< Step size is modifiable down (decrease).
};


/** \brief Convert StepType to string. */
inline
const char* toString(const StepType stepType)
{
  switch(stepType) {
    case CONSTANT_STEP_SIZE:
      return "CONSTANT_STEP_SIZE";
    case VARIABLE_STEP_SIZE:
      return "VARIABLE_STEP_SIZE";
    case UNMODIFIABLE_STEP_SIZE:
      return "UNMODIFIABLE_STEP_SIZE";
    case MODIFIABLE_STEP_SIZE:
      return "MODIFIABLE_STEP_SIZE";
    case MODIFIABLE_UP_STEP_SIZE:
      return "MODIFIABLE_UP_STEP_SIZE";
    case MODIFIABLE_DOWN_STEP_SIZE:
      return "MODIFIABLE_DOWN_STEP_SIZE";
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid IntegratorStepType!");
  }
  return 0; // Should never get here!
}


/** \brief Solution Status. */
enum SolutionStatus {
  WORKING,    ///< Solution is the working solution and not yet completed
  PASSED,     ///< Solution has passed and is archivable.
  FAILED,     ///< Solution has failed.
};


/** \brief Convert SolutionStatus to string. */
inline
const char* toString(const SolutionStatus status)
{
  switch(status) {
    case WORKING:
      return "WORKING";
    case PASSED:
      return "PASSED";
    case FAILED:
      return "FAILED";
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid SolutionStatus!");
  }
  return 0; // Should never get here!
}


} // namespace Tempus
#endif // TEMPUS_TYPES_HPP
