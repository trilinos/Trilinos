// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeStepControlStrategy_BasicBS_hpp
#define Tempus_TimeStepControlStrategy_BasicBS_hpp

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategy.hpp"
#include "Tempus_SolutionState.hpp"
#include "Tempus_SolutionStateMetaData.hpp"
#include "Tempus_SolutionHistory.hpp"
#include "Tempus_StepperState.hpp"

//Thyra
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

/** \brief StepControlStrategy class for TimeStepControl
 *
 * Section 2.2.1 / Algorithm 2.4 of A. Denner, "Experiments on
 * Temporal Variable Step BDF2 Algorithms", Masters Thesis,
 * U Wisconsin-Madison, 2014.
 *
 */
template<class Scalar>
class TimeStepControlStrategyBasicVS
  : virtual public TimeStepControlStrategy<Scalar>
{
public:

  /// Constructor
  TimeStepControlStrategyBasicVS(
    Teuchos::RCP<Teuchos::ParameterList> pList = Teuchos::null){
     this->setParameterList(pList);
  }

  /// Destructor
  virtual ~TimeStepControlStrategyBasicVS(){}

  /** \brief Determine the time step size.*/
  virtual void getNextTimeStep(
    const TimeStepControl<Scalar> tsc,
    Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory,
    Status & integratorStatus) override
  {
    using Teuchos::RCP;
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionStateMetaData<Scalar> > metaData = workingState->getMetaData();
    const Scalar errorAbs = metaData->getErrorAbs();
    const Scalar errorRel = metaData->getErrorRel();
    int order = metaData->getOrder();
    Scalar dt = metaData->getDt();
    RCP<StepperState<Scalar> > stepperState = workingState->getStepperState();
    bool printChanges = solutionHistory->getVerbLevel() !=
                        Teuchos::as<int>(Teuchos::VERB_NONE);

    RCP<Teuchos::FancyOStream> out = tsc.getOStream();
    Teuchos::OSTab ostab(out,1,"getNextTimeStep");

    auto changeDT = [] (Scalar dt_old, Scalar dt_new, std::string reason) {
      std::stringstream message;
      message <<
      "     (dt = "<<std::scientific<<std::setw(9)<<std::setprecision(3)<<dt_old
      << ", new = "<<std::scientific<<std::setw(9)<<std::setprecision(3)<<dt_new
      << ")  " << reason << std::endl;
      return message.str();
    };

    Scalar rho   = getAmplFactor();
    Scalar sigma = getReductFactor();
    Scalar eta   = computeEta(tsc, solutionHistory);

    // General rule: only increase/decrease dt once for any given reason.
    if (stepperState->stepperStatus_ == Status::FAILED) {
      if (printChanges) *out << changeDT(dt, dt*sigma,
        "Stepper failure - Decreasing dt.");
      dt *= sigma;
    }
    else { //Stepper passed
      if (eta < getMinEta()) { // increase dt
        if (printChanges) *out << changeDT(dt, dt*rho,
          "Monitoring Value (eta) is too small ("
          + std::to_string(eta) + " < " + std::to_string(getMinEta())
          + ").  Increasing dt.");
        dt *= rho;
      }
      else if (eta > getMaxEta()) { // reduce dt
        if (printChanges) *out << changeDT(dt, dt*sigma,
          "Monitoring Value (eta) is too large ("
          + std::to_string(eta) + " > " + std::to_string(getMaxEta())
          + ").  Decreasing dt.");
        dt *= sigma;
      }
      else if (errorAbs > tsc.getMaxAbsError()) { // reduce dt
        if (printChanges) *out << changeDT(dt, dt*sigma,
          "Absolute error is too large ("
          + std::to_string(errorAbs)+" > "+std::to_string(tsc.getMaxAbsError())
          + ").  Decreasing dt.");
        dt *= sigma;
      }
      else if (errorRel > tsc.getMaxRelError()) { // reduce dt
        if (printChanges) *out << changeDT(dt, dt*sigma,
          "Relative error is too large ("
          + std::to_string(errorRel)+" > "+std::to_string(tsc.getMaxRelError())
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

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
    void setParameterList(
      const Teuchos::RCP<Teuchos::ParameterList> & pList) override
    {
      if (pList == Teuchos::null) {
        // Create default parameters if null, otherwise keep current parameters.
        if (tscsPL_ == Teuchos::null) {
           tscsPL_ = Teuchos::parameterList("TimeStepControlStrategyBasicVS");
           *tscsPL_= *(this->getValidParameters());
        }
      } else {
         tscsPL_ = pList;
      }
      tscsPL_->validateParametersAndSetDefaults(*this->getValidParameters());

      TEUCHOS_TEST_FOR_EXCEPTION(getAmplFactor() <= 1.0, std::out_of_range,
      "Error - Invalid value of Amplification Factor = " << getAmplFactor()
      << "!  \n" << "Amplification Factor must be > 1.0.\n");

      TEUCHOS_TEST_FOR_EXCEPTION(getReductFactor() >= 1.0, std::out_of_range,
      "Error - Invalid value of Reduction Factor = " << getReductFactor()
      << "!  \n" << "Reduction Factor must be < 1.0.\n");

      TEUCHOS_TEST_FOR_EXCEPTION(getMinEta() > getMaxEta(), std::out_of_range,
      "Error - Invalid values of 'Minimum Value Monitoring Function' = "
      << getMinEta() << "\n and 'Maximum Value Monitoring Function' = "
      << getMaxEta() <<"! \n Mininum Value cannot be > Maximum Value! \n");

    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters() const override {
       Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

       // From (Denner, 2014), amplification factor can be at most 1.91 for
       // stability.
       pl->set<double>("Amplification Factor", 1.75, "Amplification factor");
       pl->set<double>("Reduction Factor"    , 0.5 , "Reduction factor");
       //FIXME? may need to modify default values of monitoring function
       //IKT, 1/5/17: from (Denner, 2014), it seems a reasonable choice for eta_min is 0.1*eta_max
       //Numerical tests confirm this. TODO: Change default value of eta_min to 1.0e-2?
       pl->set<double>("Minimum Value Monitoring Function" , 1.0e-6      , "Min value eta");
       pl->set<double>("Maximum Value Monitoring Function" , 1.0e-1      , "Max value eta");
       pl->set<std::string>("Name", "Basic VS");
       return pl;
    }

    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList() override {
       return tscsPL_;
    }

    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList() override {
       Teuchos::RCP<Teuchos::ParameterList> temp_plist = tscsPL_;
       tscsPL_ = Teuchos::null;
       return(temp_plist);
    }
  //@}

    virtual Scalar getAmplFactor() const
      { return tscsPL_->get<double>("Amplification Factor"); }
    virtual Scalar getReductFactor() const
      { return tscsPL_->get<double>("Reduction Factor");}
    virtual Scalar getMinEta() const
      { return tscsPL_->get<double>   ("Minimum Value Monitoring Function"); }
    virtual Scalar getMaxEta() const
      { return tscsPL_->get<double>   ("Maximum Value Monitoring Function"); }

    Scalar computeEta(const TimeStepControl<Scalar> tsc,
          const Teuchos::RCP<SolutionHistory<Scalar> > & solutionHistory)
    {
       using Teuchos::RCP;
       Scalar eta;
       const double eps = 1.0e4*std::numeric_limits<double>::epsilon();
       RCP<Teuchos::FancyOStream> out = tsc.getOStream();
       int numStates = solutionHistory->getNumStates();
       //Compute eta
       if (numStates < 3) {
          eta = getMinEta();
          return eta;
       }
       RCP<const Thyra::VectorBase<Scalar> > xOld = (*solutionHistory)[numStates-3]->getX();
       RCP<const Thyra::VectorBase<Scalar> > x = (*solutionHistory)[numStates-1]->getX();
       //IKT: uncomment the following to get some debug output
       //#define VERBOSE_DEBUG_OUTPUT
#ifdef VERBOSE_DEBUG_OUTPUT
       Teuchos::Range1D range;
       *out << "\n*** xOld ***\n";
       RTOpPack::ConstSubVectorView<Scalar> xOldv;
       xOld->acquireDetachedView(range, &xOldv);
       auto xoa = xOldv.values();
       for (auto i = 0; i < xoa.size(); ++i) *out << xoa[i] << " ";
       *out << "\n*** xOld ***\n";
       *out << "\n*** x ***\n";
       RTOpPack::ConstSubVectorView<Scalar> xv;
       x->acquireDetachedView(range, &xv);
       auto xa = xv.values();
       for (auto i = 0; i < xa.size(); ++i) *out << xa[i] << " ";
       *out << "\n*** x ***\n";
#endif
       //xDiff = x - xOld
       RCP<Thyra::VectorBase<Scalar> > xDiff = Thyra::createMember(x->space());
       Thyra::V_VmV(xDiff.ptr(), *x, *xOld);
       Scalar xDiffNorm = Thyra::norm(*xDiff);
       Scalar xOldNorm = Thyra::norm(*xOld);
       //eta = ||x^(n+1)-x^n||/(||x^n||+eps)
       eta = xDiffNorm/(xOldNorm + eps);
#ifdef VERBOSE_DEBUG_OUTPUT
       *out << "IKT xDiffNorm, xOldNorm, eta = " << xDiffNorm << ", " << xOldNorm
          << ", " << eta << "\n";
#endif
       return eta;
    }

private:
    Teuchos::RCP<Teuchos::ParameterList> tscsPL_;
};
} // namespace Tempus
#endif // Tempus_StepControlStrategy_BasicBS_hpp
