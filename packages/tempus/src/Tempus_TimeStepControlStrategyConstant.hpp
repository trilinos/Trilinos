// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControlStrategy_Constant_hpp
#define Tempus_TimeStepControlStrategy_Constant_hpp

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"


namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 */
template<class Scalar>
class TimeStepControlStrategyConstant
  : virtual public TimeStepControlStrategy<Scalar>
{
public:

  /// Constructor
  TimeStepControlStrategyConstant(){}

  /// Destructor
  virtual ~TimeStepControlStrategyConstant(){}

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(const TimeStepControl<Scalar> tsc,
    Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
    Status & integratorStatus) override
  {
     using Teuchos::RCP;
     RCP<SolutionState<Scalar> >workingState=solutionHistory->getWorkingState();
     const Scalar errorAbs = workingState->getErrorAbs();
     const Scalar errorRel = workingState->getErrorRel();
     Scalar dt = workingState->getTimeStep();

     dt = tsc.getInitTimeStep();

     RCP<Teuchos::FancyOStream> out = tsc.getOStream();
     Teuchos::OSTab ostab(out,1,"getNextTimeStep");

     // Stepper failure
     if (workingState->getSolutionStatus() == Status::FAILED) {
       *out << "Failure - Stepper failed and can not change time step size!\n"
            << "    Time step type == CONSTANT_STEP_SIZE\n" << std::endl;
       integratorStatus = FAILED;
       return;
     }

     // Absolute error failure
     if (errorAbs > tsc.getMaxAbsError()) {
       *out << "Failure - Absolute error failed and can not change time step!\n"
            << "  Time step type == CONSTANT_STEP_SIZE\n"
            << "  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="
            << tsc.getMaxAbsError() << ")" << std::endl;
       integratorStatus = FAILED;
       return;
     }

     // Relative error failure
     if (errorRel > tsc.getMaxRelError()) {
       *out << "Failure - Relative error failed and can not change time step!\n"
          << "  Time step type == CONSTANT_STEP_SIZE\n"
          << "  (errorRel ="<<errorRel<<") > (errorMaxRel ="
          << tsc.getMaxRelError() << ")" << std::endl;
       integratorStatus = FAILED;
       return;
     }

     // update dt
     workingState->setTimeStep(dt);
  }


  /// \name Overridden from Teuchos::Describable
  //@{
    std::string description() const override
    { return "Tempus::TimeStepControlStrategyConstant"; }

    void describe(Teuchos::FancyOStream          &out,
                  const Teuchos::EVerbosityLevel verbLevel) const override
    {
      Teuchos::OSTab ostab(out,2,"describe");
      out << description() << std::endl;
    }
  //@}

};


} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_Constant_hpp
