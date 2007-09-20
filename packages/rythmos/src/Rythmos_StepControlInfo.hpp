
#ifndef RYTHMOS_STEP_CONTROL_INFO_HPP
#define RYTHMOS_STEP_CONTROL_INFO_HPP

#include "Rythmos_StepperSupportTypes.hpp"

namespace Rythmos {


/** \brief Simple strict to aggregate integration/stepper control information.
 */
template<class Scalar>
struct StepControlInfo {
  /** \brief The size of the time step. */
  Scalar stepSize;
  /** \brief The type of time step. */
  StepSizeType stepType;
  /** \brief True if step size is limited by a breakpoint. */
  bool limitedByBreakPoint;
  /** \brief True if the time integrator should restart when passing
   * over the breakpoint. */
  EBreakPointType breakPointType;
  /** \brief Initalize to invalid state. */
  StepControlInfo()
    :stepSize(-1.0), stepType(STEP_TYPE_FIXED),
     limitedByBreakPoint(false),
     breakPointType(BREAK_POINT_TYPE_SOFT)
    {}
};


/** \brief . */
template<class Scalar>
std::ostream& operator<<(
  std::ostream &out, const StepControlInfo<Scalar> &stepCtrlInfo
  )
{
  using std::endl;
  out
    << "stepType = " << toString(stepCtrlInfo.stepType) << endl
    << "stepSize = " << stepCtrlInfo.stepSize << endl
    << "limitedByBreakPoint = " << stepCtrlInfo.limitedByBreakPoint << endl
    << "breakPointType = " << toString(stepCtrlInfo.breakPointType) << endl;
  return out;
}


/** \brief Return a step control info object for a step actually taken.
 *
 * All integrator subclass implementations should use this function to update
 * the step control info for a step actually taken.
 *
 * \relates StepControlInfo
 */
template<class Scalar>
StepControlInfo<Scalar>
stepCtrlInfoTaken( 
  const StepControlInfo<Scalar> &trialStepCtrlInfo,
  const Scalar &stepSizeTaken
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  const Scalar zero = ST::zero();
  StepControlInfo<Scalar> stepCtrlInfo = trialStepCtrlInfo;
  stepCtrlInfo.stepSize = stepSizeTaken;
  if ( trialStepCtrlInfo.stepSize > zero && stepSizeTaken > zero ) {
    if (stepSizeTaken < trialStepCtrlInfo.stepSize) {
      stepCtrlInfo.limitedByBreakPoint = false;
    }
  }
  else {
    stepCtrlInfo.limitedByBreakPoint = false;
  }
  return stepCtrlInfo;
}


// 2007/09/14: rabartl: ToDo: Above: Move this function into
// Rythmos_IntegratorBaseHelpers.hpp once created!


} // namespace Rythmos


#endif // RYTHMOS_STEP_CONTROL_INFO_HPP
