// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepControlStrategy_VSBDF2_hpp
#define Tempus_StepControlStrategy_VSBDF2_hpp

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_StepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"


namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 * Section 2.2.1 / Algorithm 2.4 of A. Denner, "Experiments on
 * Temporal Variable Step BDF2 Algorithms", Masters Thesis,
 * U Wisconsin-Madison, 2014.
 *
 */
template<class Scalar>
class StepControlStrategyVSBDF2 : virtual public StepControlStrategy<Scalar>
{
public:

  /// Constructor
  StepControlStrategyVSBDF2(){}

  /// Destructor
  virtual ~StepControlStrategyVSBDF2(){}

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(TimeStepControl<Scalar> tsc, Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
        Status & integratorStatus) override {

     Teuchos::RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
     Teuchos::RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
     const Scalar errorAbs = metaData->getErrorAbs();
     const Scalar errorRel = metaData->getErrorRel();
     int order = metaData->getOrder();
     Scalar dt = metaData->getDt();
     Teuchos::RCP<StepperState<Scalar> > stepperState = workingState->getStepperState();
     bool printChanges = solutionHistory->getVerbLevel() !=
                        Teuchos::as<int>(Teuchos::VERB_NONE);

     Teuchos::RCP<Teuchos::FancyOStream> out = tsc.getOStream();
    Teuchos::OSTab ostab(out,1,"getNextTimeStep");

    //TODO: prevent redefining this
    auto changeDT = [] (Scalar dt_old, Scalar dt_new, std::string reason) {
       std::stringstream message;
       message <<
          "     (dt = "<<std::scientific<<std::setw(9)<<std::setprecision(3)<<dt_old
          << ", new = "<<std::scientific<<std::setw(9)<<std::setprecision(3)<<dt_new
          << ")  " << reason << std::endl;
       return message.str();
    };

      Scalar rho   = tsc.getAmplFactor();
      Scalar sigma = tsc.getReductFactor();
      Scalar eta   = tsc.computeEta(solutionHistory);

      // General rule: only increase/decrease dt once for any given reason.
      if (stepperState->stepperStatus_ == Status::FAILED) {
        if (printChanges) *out << changeDT(dt, dt*sigma,
          "Stepper failure - Decreasing dt.");
        dt *= sigma;
      }
      else { //Stepper passed
        if (eta < tsc.getMinEta()) { // increase dt
          if (printChanges) *out << changeDT(dt, dt*rho,
            "Monitoring Value (eta) is too small ("
            + std::to_string(eta) + " < " + std::to_string(tsc.getMinEta())
            + ").  Increasing dt.");
          dt *= rho;
        }
        else if (eta > tsc.getMaxEta()) { // reduce dt
          if (printChanges) *out << changeDT(dt, dt*sigma,
            "Monitoring Value (eta) is too large ("
            + std::to_string(eta) + " > " + std::to_string(tsc.getMaxEta())
            + ").  Decreasing dt.");
          dt *= sigma;
        }
        else if (errorAbs > tsc.getMaxAbsError()) { // reduce dt
          if (printChanges) *out << changeDT(dt, dt*sigma,
            "Absolute error is too large ("
            + std::to_string(errorAbs) +" > "+ std::to_string(tsc.getMaxAbsError())
            + ").  Decreasing dt.");
          dt *= sigma;
        }
        else if (errorRel > tsc.getMaxRelError()) { // reduce dt
          if (printChanges) *out << changeDT(dt, dt*sigma,
            "Relative error is too large ("
            + std::to_string(errorRel) +" > "+ std::to_string(tsc.getMaxRelError())
            + ").  Decreasing dt.");
          dt *= sigma;
        }
        else if (order < tsc.getMinOrder()) { // order too low, increase dt
          if (printChanges) *out << changeDT(dt, dt*rho,
            "Order is too small ("
            + std::to_string(order) + " < " + std::to_string(tsc.getMinOrder())
            + ").  Increasing dt.");
          dt *= rho;
        }
        else if (order > tsc.getMaxOrder()) { // order too high, reduce dt
          if (printChanges) *out << changeDT(dt, dt*sigma,
            "Order is too large ("
            + std::to_string(order) + " > " + std::to_string(tsc.getMaxOrder())
            + ").  Decreasing dt.");
          dt *= sigma;
        }
      }

      if (dt < tsc.getMinTimeStep()) { // decreased below minimum dt
        if (printChanges) *out << changeDT(dt, tsc.getMinTimeStep(),
          "dt is too small.  Resetting to minimum dt.");
        dt = tsc.getMinTimeStep();
      }
      if (dt > tsc.getMaxTimeStep()) { // increased above maximum dt
        if (printChanges) *out << changeDT(dt, tsc.getMaxTimeStep(),
          "dt is too large.  Resetting to maximum dt.");
        dt = tsc.getMaxTimeStep();
      }

      metaData->setOrder(order);
      metaData->setDt(dt);
    }

};
} // namespace Tempus
#endif // Tempus_StepControlStrategy_VSBDF2_hpp
