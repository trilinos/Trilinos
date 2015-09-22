#include "Pike_Rxn_ModelEvaluator_All.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <cmath>

namespace pike_test {

  RxnAll::RxnAll(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd) :
    mpd_(mpd),
    x_(3),
    xOld_(3),
    xAnalytic_(3),
    Y_(3),
    f_(3),
    p_k1_(1.0),
    p_k2_(2.0),
    p_CA0_(1.0),
    p_CB0_(0.0),
    p_CC0_(0.0),
    timeStepNumber_(1),
    isLocallyConverged_(false)
  {
    this->reset();
  }

  std::string RxnAll::name() const
  { return "EqAll"; }
  
  void RxnAll::solve()
  {
    tentativeTime_ = this->getCurrentTime() + this->getCurrentTimeStepSize();

    // Explicit RK4

    // s=1
    double t = this->getCurrentTime();
    for (int i=0; i<3; ++i)
      Y_[i] = xOld_[i]; // Y1
    this->evaluateF(t,Y_,f_);
    for (int i=0; i<3; ++i)
      x_[i] = xOld_[i] + this->getCurrentTimeStepSize() / 6.0 * f_[i];

    // s=2
    for (int i=0; i<3; ++i)
      Y_[i] = xOld_[i] + 0.5 * this->getCurrentTimeStepSize() * f_[i]; //Y2
    t = this->getCurrentTime() + 0.5 * this->getCurrentTimeStepSize();
    this->evaluateF(t,Y_,f_);
    for (int i=0; i<3; ++i)
      x_[i] += this->getCurrentTimeStepSize() / 3.0 * f_[i];
    
    // s=3
    for (int i=0; i<3; ++i)
      Y_[i] = xOld_[i] + 0.5 * this->getCurrentTimeStepSize() * f_[i]; // Y3
    t = this->getCurrentTime() + 0.5 * this->getCurrentTimeStepSize();
    this->evaluateF(t,Y_,f_);
    for (int i=0; i <3; ++i)
      x_[i] += this->getCurrentTimeStepSize() / 3.0 * f_[i];
    
    // s=4
    for (int i=0; i <3; ++i)
      Y_[i] = xOld_[i] + this->getCurrentTimeStepSize() * f_[i]; // Y4
    t = this->getCurrentTime() + this->getCurrentTimeStepSize();
    this->evaluateF(t,Y_,f_);
    for (int i=0; i <3; ++i)
      x_[i] += this->getCurrentTimeStepSize() / 6.0 * f_[i];


    //std::cout << "time=" << t << ", x=" << x_[0] << "," << x_[1] << "," <<  x_[2] << std::endl;

    // Teuchos::broadcast(*mpd_->getGlobalComm(), 0, const_cast<double*>(&x1_));

    isLocallyConverged_ = true;

    solvedTentativeStep_ = true;
  }
  
  bool RxnAll::isLocallyConverged() const
  {
    return isLocallyConverged_;
  }

  bool RxnAll::supportsParameter(const std::string& pName) const
  {
    return false;
  }

  int RxnAll::getNumberOfParameters() const
  {
    return 0;
  }

  std::string RxnAll::getParameterName(const int l) const
  {
    TEUCHOS_ASSERT(false);
    return "No supported parameters!";
  }

  int RxnAll::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_ASSERT(false);
    return -1;
  }

  void RxnAll::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_ASSERT(false);
  }

  Teuchos::ArrayView<const double> RxnAll::getResponse(const int i) const
  {
    TEUCHOS_ASSERT(false);
    return Teuchos::ArrayView<const double>(&(x_[0]),1);
  }
  
  int RxnAll::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_ASSERT(false);
    return 0;
  }
  
  std::string RxnAll::getResponseName(const int i) const
  {
    TEUCHOS_ASSERT(false);
    return "No responses supported!";
  }

  bool RxnAll::supportsResponse(const std::string& rName) const
  {
    return false;
  }
  
  int RxnAll::getNumberOfResponses() const
  {
    return 0;
  }

  bool RxnAll::isTransient() const
  {
    return true;
  }

  double RxnAll::getCurrentTime() const
  {
    return currentTime_;
  }

  double RxnAll::getTentativeTime() const
  {
    return tentativeTime_;
  }

  bool RxnAll::solvedTentativeStep() const
  {
    return solvedTentativeStep_;
  }

  double RxnAll::getCurrentTimeStepSize() const
  {
    return currentTimeStepSize_;
  }
  
  double RxnAll::getDesiredTimeStepSize() const
  {
    return 1.0;
  }

  double RxnAll::getMaxTimeStepSize() const
  {
    return 10.0;
  }
  
  void RxnAll::setNextTimeStepSize(const double& dt)
  {
    currentTimeStepSize_ = dt;
  }

  void RxnAll::acceptTimeStep() 
  {
    ++timeStepNumber_;
    currentTime_ = tentativeTime_;
    solvedTentativeStep_ = false;

    for (int i=0; i<3; ++i)
      xOld_[i] = x_[i];
  }

  void RxnAll::evaluateF(const double& ,
			 const std::vector<double>& x, 
			 std::vector<double>& f)
  {
    // Simple parallel rxn for A->B at rate k1 and A->C at rate k2
    f[0] = -p_k1_ * x[0] - p_k2_ * x[0];
    f[1] = +p_k1_ * x[0];
    f[2] = +p_k2_ * x[0];
  }

  double RxnAll::evaluateError()
  {
    xAnalytic_[0] = p_CA0_ * std::exp(- p_k1_ * currentTime_ - p_k2_ * currentTime_);
    xAnalytic_[1] = p_CB0_ + p_CA0_ * p_k1_ / (p_k1_ + p_k2_) * (1.0 - std::exp(-(p_k1_+p_k2_)*currentTime_));
    xAnalytic_[2] = p_CC0_ + p_CA0_ * p_k2_ / (p_k1_ + p_k2_) * (1.0 - std::exp(-(p_k1_+p_k2_)*currentTime_));

    //std::cout << "t=" << currentTime_ << ", x=" << x_[0] << "," << x_[1] << "," <<  x_[2] << std::endl;

    double error = 0.0;
    for (int i=0; i<3; ++i)
      error += (x_[i] - xAnalytic_[i]) * (x_[i] - xAnalytic_[i]);
    return std::sqrt(error);
  }

  void RxnAll::reset()
  {
    timeStepNumber_ = 0;
    solvedTentativeStep_ = false;
    currentTime_ = 0.0;
    tentativeTime_ = 0.0;
    currentTimeStepSize_ = 1.0;
    isLocallyConverged_ = false;

    xOld_[0] = p_CA0_;
    x_[0] = p_CA0_;
    xOld_[1] = p_CB0_;
    x_[1] = p_CB0_;
    xOld_[2] = p_CC0_;
    x_[2] = p_CC0_;
  }

  // non-member ctor
  Teuchos::RCP<pike_test::RxnAll> 
  rxnAll(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd)
  { return Teuchos::rcp(new pike_test::RxnAll(mpd)); }
  
}
