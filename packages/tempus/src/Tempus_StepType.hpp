#ifndef TEMPUS_STEPTYPE_HPP
#define TEMPUS_STEPTYPE_HPP

namespace tempus {


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
  PASSING,     ///< Constant integrator step size
  FAILED,     ///< Variable integrator step size
};


/** \brief Convert StepType to string. */
inline
const char* toString(const SolutionStatus status)
{
  switch(status) {
    case PASSING:
      return "PASSING";
    case FAILED:
      return "FAILED";
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid SolutionStatus!");
  }
  return 0; // Should never get here!
}


} // namespace tempus
#ifndef TEMPUS_STEPTYPE_HPP
