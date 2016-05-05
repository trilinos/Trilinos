#ifndef TEMPUS_STEPTYPE_HPP
#define TEMPUS_STEPTYPE_HPP

namespace Tempus {


/** \brief Step type. */
enum StepType {
  CONSTANT_STEP_SIZE,     ///< Constant integrator step size
  VARIABLE_STEP_SIZE,     ///< Variable integrator step size
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
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid StepType!");
  }
  return 0; // Should never get here!
}


/** \brief Solution Status. */
enum SolutionStatus {
  WORKING,    ///< Solution is the working solution and not yet completed
  PASSED,     ///< Solution has passed and is archivable.
  FAILED,     ///< Solution has failed.
};


/** \brief Convert StepType to string. */
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
#endif // TEMPUS_STEPTYPE_HPP
