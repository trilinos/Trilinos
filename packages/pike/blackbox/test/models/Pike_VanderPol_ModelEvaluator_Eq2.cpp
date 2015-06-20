#include "Pike_VanderPol_ModelEvaluator_Eq2.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <cmath>

namespace pike_test {

  VanderPolME2::VanderPolME2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd) :
    mpd_(mpd),
    eq2_(mpd->getApplicationIndex("Eq2")),
    x2_(0.0),
    x2old_(x2_),
    p_epsilon_(1.0e-3),
    p_x1_(2.0),
    isLocallyConverged_(true)
  {
    
  }

  std::string VanderPolME2::name() const
  { return "Eq2"; }
  
  double VanderPolME2::computeF(double x2)
  {
    return (x2_ - x2old_) / currentTimeStepSize_ + p_epsilon_ * (1.0 - p_x1_ * p_x1_ * x2_ - p_x1_);
  }

  double VanderPolME2::computeJ(double x2)
  {
    return 1.0 / currentTimeStepSize_ - p_epsilon_ * p_x1_ * p_x1_;
  }

  void VanderPolME2::solve()
  {
    if (mpd_->appExistsOnProcess(eq2_)) {

      // take a tentative time step
      TEUCHOS_ASSERT(true);
      tentativeTime_ = currentTime_ + currentTimeStepSize_;
      
      // Solve for x2 at new time using Newton's method with backward
      // Euler.
      //const double& t = tentativeTime_;
      const double relTol = 1.0e-5;
      const double absTol = 1.0e-8;
      bool converged = false;
      
      const int maxIters = 20;
      int iter = 0;
      double f = computeF(x2_);
      *mpd_->getApplicationOStream(eq2_) << "Solving for x2_ in me2, x1 = " << p_x1_ << std::endl;
      *mpd_->getApplicationOStream(eq2_) << iter << ": x2=" << x2_ << ", ||f||=" << std::abs(f) << std::endl;
      while (!converged && (iter < maxIters) ) {
	
	double dx = -f / computeJ(x2_);
	x2_ = dx + x2_;
	f = computeF(x2_);
	++iter;
	
	if ( (std::abs(dx) < relTol) && (std::abs(f) < absTol) )
	  converged = true;
	
	*mpd_->getApplicationOStream(eq2_) << iter << ": x2=" << x2_ << ", ||f||=" << std::abs(f) << std::endl;
      }

      if (converged)
	isLocallyConverged_ = true;
      else
	isLocallyConverged_ = false; 
      
      *mpd_->getApplicationOStream(eq2_) << "me2::solve(), x2_ = " << x2_ << std::endl; 
    }

    Teuchos::broadcast(*mpd_->getGlobalComm(), 1, const_cast<double*>(&x2_));

    // Send convergence status from proc 1 where eq 2 lives to all
    // procs.  Use int for mpi commmunication.
    int zeroMeansConverged = isLocallyConverged_ ? 0 : 1;
    Teuchos::broadcast(*mpd_->getGlobalComm(), 1, Teuchos::ptrFromRef(zeroMeansConverged));
    isLocallyConverged_ = (zeroMeansConverged == 0);

    sovledTentativeStep_ = true;
  }
  
  bool VanderPolME2::isLocallyConverged() const
  {
    return isLocallyConverged_;
  }

  bool VanderPolME2::supportsParameter(const std::string& pName) const
  {
    if ( (pName == "epsilon") || (pName == "x1") )
      return true;
    
    return false;
  }

  int VanderPolME2::getNumberOfParameters() const
  {
    return 2;
  }

  std::string VanderPolME2::getParameterName(const int l) const
  {
    TEUCHOS_ASSERT(l >= 0);
    TEUCHOS_ASSERT(l < 2);
    if (l == 0)
      return "epsilon";
    
    return "x1";
  }

  int VanderPolME2::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_ASSERT( (pName == "epsilon") || (pName == "x1") );
    
    if (pName == "epsilon")
      return 0;

    return 1;
  }

  void VanderPolME2::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_ASSERT(l >= 0);
    TEUCHOS_ASSERT(l < 2);
    if (l == 0)
      p_epsilon_ = p[0];
    
    p_x1_ = p[0];
  }

  Teuchos::ArrayView<const double> VanderPolME2::getResponse(const int i) const
  {
    TEUCHOS_ASSERT(i == 0);
    return Teuchos::ArrayView<const double>(&x2_,1);
  }
  
  int VanderPolME2::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_ASSERT(rName == "x2");
    return 0;
  }
  
  std::string VanderPolME2::getResponseName(const int i) const
  {
    TEUCHOS_ASSERT(i == 0);
    return "x2";
  }

  bool VanderPolME2::supportsResponse(const std::string& rName) const
  {
    if (rName == "x2")
      return true;
    
    return false;
  }
  
  int VanderPolME2::getNumberOfResponses() const
  {
    return 1;
  }

  bool VanderPolME2::isTransient() const
  {
    return true;
  }

  double VanderPolME2::getCurrentTime() const
  {
    return currentTime_;
  }

  double VanderPolME2::getTentativeTime() const
  {
    return tentativeTime_;
  }

  bool VanderPolME2::solvedTentativeStep() const
  {
    return sovledTentativeStep_;
  }

  double VanderPolME2::getCurrentTimeStepSize() const
  {
    return currentTimeStepSize_;
  }
  
  double VanderPolME2::getDesiredTimeStepSize() const
  {
    return 1.0;
  }

  double VanderPolME2::getMaxTimeStepSize() const
  {
    return 10.0;
  }
  
  void VanderPolME2::setNextTimeStepSize(const double& dt)
  {
    currentTimeStepSize_ = dt;
  }

  void VanderPolME2::acceptTimeStep() 
  {
    x2old_ = x2_;
    currentTime_ = tentativeTime_;
    sovledTentativeStep_ = false;
  }
  
  // non-member ctor
  Teuchos::RCP<pike_test::VanderPolME2> 
  vanderPolME2(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd)
  { return Teuchos::rcp(new pike_test::VanderPolME2(mpd)); }
  
}
