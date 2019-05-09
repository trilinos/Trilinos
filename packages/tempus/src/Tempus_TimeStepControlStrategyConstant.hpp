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
#include "Tempus_SolutionStateMetaData.hpp"
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
     RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
     const Scalar errorAbs = metaData->getErrorAbs();
     const Scalar errorRel = metaData->getErrorRel();
     int order = metaData->getOrder();
     Scalar dt = metaData->getDt();
     bool printChanges = solutionHistory->getVerbLevel() !=
        Teuchos::as<int>(Teuchos::VERB_NONE);

     dt = tsc.getInitTimeStep();

     RCP<Teuchos::FancyOStream> out = tsc.getOStream();
     Teuchos::OSTab ostab(out,1,"getNextTimeStep");

     auto changeOrder = [] (int order_old, int order_new, std::string reason) {
        std::stringstream message;
        message << "     (order = " << std::setw(2) << order_old
           <<       ", new = " << std::setw(2) << order_new
           << ")  " << reason << std::endl;
        return message.str();
     };

     // Stepper failure
     if (workingState->getSolutionStatus() == Status::FAILED) {
        if (order+1 <= tsc.getMaxOrder()) {
           if (printChanges) *out << changeOrder(order, order+1,
                 "Stepper failure, increasing order.");
           order++;
        } else {
           *out << "Failure - Stepper failed and can not change time step size "
              << "or order!\n"
              << "    Time step type == CONSTANT_STEP_SIZE\n"
              << "    order = " << order << std::endl;
           integratorStatus = FAILED;
           return;
        }
     }

     // Absolute error failure
     if (errorAbs > tsc.getMaxAbsError()) {
        if (order+1 <= tsc.getMaxOrder()) {
           if (printChanges) *out << changeOrder(order, order+1,
                 "Absolute error is too large.  Increasing order.");
           order++;
        } else {
           *out
              << "Failure - Absolute error failed and can not change time step "
              << "size or order!\n"
              << "  Time step type == CONSTANT_STEP_SIZE\n"
              << "  order = " << order
              << "  (errorAbs ="<<errorAbs<<") > (errorMaxAbs ="
              << tsc.getMaxAbsError()<<")"
              << std::endl;
           integratorStatus = FAILED;
           return;
        }
     }

     // Relative error failure
     if (errorRel > tsc.getMaxRelError()) {
        if (order+1 <= tsc.getMaxOrder()) {
           if (printChanges) *out << changeOrder(order, order+1,
                 "Relative error is too large.  Increasing order.");
           order++;
        } else {
           *out
              << "Failure - Relative error failed and can not change time step "
              << "size or order!\n"
              << "  Time step type == CONSTANT_STEP_SIZE\n"
              << "  order = " << order
              << "  (errorRel ="<<errorRel<<") > (errorMaxRel ="
              << tsc.getMaxRelError()<<")"
              << std::endl;
           integratorStatus = FAILED;
           return;
        }
     }

     // Consistency checks
     TEUCHOS_TEST_FOR_EXCEPTION(
       (order < tsc.getMinOrder() || order > tsc.getMaxOrder()),
       std::out_of_range,
       "Error - Solution order is out of range and can not change "
       "time step size!\n"
       "    Time step type == CONSTANT_STEP_SIZE\n"
       "    [order_min, order_max] = [" <<tsc.getMinOrder()<< ", "
       <<tsc.getMaxOrder()<< "]\n"
       "    order = " << order << "\n");

     // update order and dt
     metaData->setOrder(order);
     metaData->setDt(dt);
  }

};
} // namespace Tempus
#endif // Tempus_TimeStepControlStrategy_Constant_hpp
