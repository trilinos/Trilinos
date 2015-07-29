#include "Pike_Rxn_ModelEvaluator_SingleEqBase.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <cmath>

namespace pike_test {

  RxnSingleEqBase::RxnSingleEqBase(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd) :
    mpd_(mpd),
    x_(1),
    xOld_(1),
    xAnalytic_(1),
    Y_(1),
    f_(1),
    p_k1_(1.0),
    p_k2_(2.0),
    p_CA0_(1.0),
    p_CB0_(0.0),
    p_CC0_(0.0),
    timeStepNumber_(1),
    isLocallyConverged_(false)
  {
    //this->reset();
  }

  std::string RxnSingleEqBase::name() const
  { return "EqAll"; }
  
  void RxnSingleEqBase::solve()
  {
    tentativeTime_ = this->getCurrentTime() + this->getCurrentTimeStepSize();

    // Explicit RK4
    const int numEqns = 1;

    // s=1
    double t = this->getCurrentTime();
    for (int i=0; i<numEqns; ++i)
      Y_[i] = xOld_[i]; // Y1
    this->evaluateF(t,Y_,f_);
    for (int i=0; i<numEqns; ++i)
      x_[i] = xOld_[i] + this->getCurrentTimeStepSize() / 6.0 * f_[i];

    // s=2
    for (int i=0; i<numEqns; ++i)
      Y_[i] = xOld_[i] + 0.5 * this->getCurrentTimeStepSize() * f_[i]; //Y2
    t = this->getCurrentTime() + 0.5 * this->getCurrentTimeStepSize();
    this->evaluateF(t,Y_,f_);
    for (int i=0; i<numEqns; ++i)
      x_[i] += this->getCurrentTimeStepSize() / 3.0 * f_[i];
    
    // s=3
    for (int i=0; i<numEqns; ++i)
      Y_[i] = xOld_[i] + 0.5 * this->getCurrentTimeStepSize() * f_[i]; // Y3
    t = this->getCurrentTime() + 0.5 * this->getCurrentTimeStepSize();
    this->evaluateF(t,Y_,f_);
    for (int i=0; i <numEqns; ++i)
      x_[i] += this->getCurrentTimeStepSize() / 3.0 * f_[i];
    
    // s=4
    for (int i=0; i <numEqns; ++i)
      Y_[i] = xOld_[i] + this->getCurrentTimeStepSize() * f_[i]; // Y4
    t = this->getCurrentTime() + this->getCurrentTimeStepSize();
    this->evaluateF(t,Y_,f_);
    for (int i=0; i <numEqns; ++i)
      x_[i] += this->getCurrentTimeStepSize() / 6.0 * f_[i];


    //std::cout << "time=" << t << ", x=" << x_[0] << "," << x_[1] << "," <<  x_[2] << std::endl;

    // Teuchos::broadcast(*mpd_->getGlobalComm(), 0, const_cast<double*>(&x1_));

    isLocallyConverged_ = true;

    solvedTentativeStep_ = true;
  }
  
  bool RxnSingleEqBase::isLocallyConverged() const
  {
    return isLocallyConverged_;
  }

  bool RxnSingleEqBase::isTransient() const
  {
    return true;
  }

  double RxnSingleEqBase::getCurrentTime() const
  {
    return currentTime_;
  }

  double RxnSingleEqBase::getTentativeTime() const
  {
    return tentativeTime_;
  }

  bool RxnSingleEqBase::solvedTentativeStep() const
  {
    return solvedTentativeStep_;
  }

  double RxnSingleEqBase::getCurrentTimeStepSize() const
  {
    return currentTimeStepSize_;
  }
  
  double RxnSingleEqBase::getDesiredTimeStepSize() const
  {
    return 1.0;
  }

  double RxnSingleEqBase::getMaxTimeStepSize() const
  {
    return 10.0;
  }
  
  void RxnSingleEqBase::setNextTimeStepSize(const double& dt)
  {
    currentTimeStepSize_ = dt;
  }

  void RxnSingleEqBase::acceptTimeStep() 
  {
    ++timeStepNumber_;
    currentTime_ = tentativeTime_;
    solvedTentativeStep_ = false;

    for (int i=0; i<1; ++i)
      xOld_[i] = x_[i];
  }
  
}
