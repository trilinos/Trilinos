#include "Pike_VanderPol_ModelEvaluator_Eq1.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <cmath>

namespace pike_test {

  VanderPolME1::VanderPolME1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd) :
    mpd_(mpd),
    eq1_(mpd->getApplicationIndex("Eq1")),
    x1_(2.0),
    p_epsilon_(1.0e-3),
    p_x2_(0.0),
    timeStepNumber_(1),
    isLocallyConverged_(false),
    unitTestSolveFailureOnTimeStepNumber_(-1),
    unitTestFailureTriggered_(false)
  {
    
  }

  std::string VanderPolME1::name() const
  { return "Eq1"; }
  
  void VanderPolME1::solve()
  {
    // This block is for unit testing
    if ( (timeStepNumber_ == unitTestSolveFailureOnTimeStepNumber_) &&
     	 (!unitTestFailureTriggered_) ) {
      isLocallyConverged_ = false;
      sovledTentativeStep_ = true;
      unitTestFailureTriggered_ = true;
      *mpd_->getApplicationOStream(eq1_) << "UNIT TEST TRIGGERING FAILURE IN ME 1" << std::endl;
      return;
    }

    if (mpd_->appExistsOnProcess(eq1_)) {

      // take a tentative time step
      TEUCHOS_ASSERT(true);
      tentativeTime_ = currentTime_ + currentTimeStepSize_;
      
      // Solve for x1 at new time.  Simple analytic solution.
      x1_ =  p_x2_ * currentTimeStepSize_  + x1_;

      *mpd_->getApplicationOStream(eq1_) << "me1::solve(), currentTimeStepSize_= " << currentTimeStepSize_ 
					 << ", p_x2_ = " << p_x2_ << std::endl; 
      *mpd_->getApplicationOStream(eq1_) << "me1::solve(), x1_ = " << x1_ << std::endl; 
    }

    Teuchos::broadcast(*mpd_->getGlobalComm(), 0, const_cast<double*>(&x1_));

    isLocallyConverged_ = true;

    sovledTentativeStep_ = true;
  }
  
  bool VanderPolME1::isLocallyConverged() const
  {
    return isLocallyConverged_;
  }

  bool VanderPolME1::supportsParameter(const std::string& pName) const
  {
    if ( (pName == "epsilon") || (pName == "x2") )
      return true;
    
    return false;
  }

  int VanderPolME1::getNumberOfParameters() const
  {
    return 2;
  }

  std::string VanderPolME1::getParameterName(const int l) const
  {
    TEUCHOS_ASSERT(l >= 0);
    TEUCHOS_ASSERT(l < 2);
    if (l == 0)
      return "epsilon";
    
    return "x2";
  }

  int VanderPolME1::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_ASSERT( (pName == "epsilon") || (pName == "x2") );
    
    if (pName == "epsilon")
      return 0;

    return 1;
  }

  void VanderPolME1::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_ASSERT(l >= 0);
    TEUCHOS_ASSERT(l < 2);
    if (l == 0)
      p_epsilon_ = p[0];
    
    p_x2_ = p[0];
  }

  Teuchos::ArrayView<const double> VanderPolME1::getResponse(const int i) const
  {
    TEUCHOS_ASSERT(i == 0);
    return Teuchos::ArrayView<const double>(&x1_,1);
  }
  
  int VanderPolME1::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_ASSERT(rName == "x1");
    return 0;
  }
  
  std::string VanderPolME1::getResponseName(const int i) const
  {
    TEUCHOS_ASSERT(i == 0);
    return "x1";
  }

  bool VanderPolME1::supportsResponse(const std::string& rName) const
  {
    if (rName == "x1")
      return true;
    
    return false;
  }
  
  int VanderPolME1::getNumberOfResponses() const
  {
    return 1;
  }

  bool VanderPolME1::isTransient() const
  {
    return true;
  }

  double VanderPolME1::getCurrentTime() const
  {
    return currentTime_;
  }

  double VanderPolME1::getTentativeTime() const
  {
    return tentativeTime_;
  }

  bool VanderPolME1::solvedTentativeStep() const
  {
    return sovledTentativeStep_;
  }

  double VanderPolME1::getCurrentTimeStepSize() const
  {
    return currentTimeStepSize_;
  }
  
  double VanderPolME1::getDesiredTimeStepSize() const
  {
    return 1.0;
  }

  double VanderPolME1::getMaxTimeStepSize() const
  {
    return 10.0;
  }
  
  void VanderPolME1::setNextTimeStepSize(const double& dt)
  {
    currentTimeStepSize_ = dt;
  }

  void VanderPolME1::acceptTimeStep() 
  {
    ++timeStepNumber_;
    currentTime_ = tentativeTime_;
    sovledTentativeStep_ = false;
  }
  
  void VanderPolME1::setUnitTestFailureForTimeStep(const int& timeStepNumber)
  {
    unitTestSolveFailureOnTimeStepNumber_ = timeStepNumber;
  }

  // non-member ctor
  Teuchos::RCP<pike_test::VanderPolME1> 
  vanderPolME1(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd)
  { return Teuchos::rcp(new pike_test::VanderPolME1(mpd)); }
  
}
