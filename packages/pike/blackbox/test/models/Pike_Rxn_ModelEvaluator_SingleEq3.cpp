#include "Pike_Rxn_ModelEvaluator_SingleEq3.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <cmath>

namespace pike_test {

  RxnSingleEq3::RxnSingleEq3(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd) :
    RxnSingleEqBase(mpd)
  {
    this->reset();
  }

  std::string RxnSingleEq3::name() const
  {
    return "Eq3";
  }

  bool RxnSingleEq3::supportsParameter(const std::string& pName) const
  {
    if ( (pName == "CA") || (pName == "CB") )
      return true;    
    return false;
  }

  int RxnSingleEq3::getNumberOfParameters() const
  {
    return 2;
  }

  std::string RxnSingleEq3::getParameterName(const int l) const
  {
    TEUCHOS_ASSERT(l >= 0);
    TEUCHOS_ASSERT(l < 2);
    if (l == 0)
      return "CA";
    return "CB";
  }

  int RxnSingleEq3::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_ASSERT( (pName == "CA") || (pName == "CB") );
    
    if (pName == "CA")
      return 0;

    return 1;
  }

  void RxnSingleEq3::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_ASSERT(l >= 0);
    TEUCHOS_ASSERT(l < 2);
    if (l == 0)
      p_CA_ = p[0];
    
    p_CB_ = p[0];
  }

  Teuchos::ArrayView<const double> RxnSingleEq3::getResponse(const int i) const
  {
    TEUCHOS_ASSERT(i == 0);
    return Teuchos::ArrayView<const double>(&(x_[0]),1);
  }
  
  int RxnSingleEq3::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_ASSERT(rName == "CC");
    return 0;
  }
  
  std::string RxnSingleEq3::getResponseName(const int i) const
  {
    TEUCHOS_ASSERT(i == 0);
    return "CC";
  }

  bool RxnSingleEq3::supportsResponse(const std::string& rName) const
  {
    if (rName == "CC")
      return true;
    
    return false;
  }
  
  int RxnSingleEq3::getNumberOfResponses() const
  {
    return 1;
  }


  void RxnSingleEq3::evaluateF(const double& ,
			 const std::vector<double>& x, 
			 std::vector<double>& f)
  {
    // Simple parallel rxn for A->B at rate k1 and A->C at rate k2
    //f[0] = -p_k1_ * x[0] - p_k2_ * x[0];
    //f[1] = +p_k1_ * p_CA_;
    f[2] = +p_k2_ * p_CA_;
  }

  double RxnSingleEq3::evaluateError()
  {
    //xAnalytic_[0] = p_CA0_ * std::exp(- p_k1_ * currentTime_ - p_k2_ * currentTime_);
    //xAnalytic_[1] = p_CB0_ + p_CA0_ * p_k1_ / (p_k1_ + p_k2_) * (1.0 - std::exp(-(p_k1_+p_k2_)*currentTime_));
    xAnalytic_[2] = p_CC0_ + p_CA0_ * p_k2_ / (p_k1_ + p_k2_) * (1.0 - std::exp(-(p_k1_+p_k2_)*currentTime_));

    //std::cout << "t=" << currentTime_ << ", x=" << x_[0] << "," << x_[1] << "," <<  x_[2] << std::endl;

    double error = 0.0;
    for (int i=0; i<1; ++i)
      error += (x_[i] - xAnalytic_[i]) * (x_[i] - xAnalytic_[i]);
    return std::sqrt(error);
  }

  void RxnSingleEq3::reset()
  {
    timeStepNumber_ = 0;
    solvedTentativeStep_ = false;
    currentTime_ = 0.0;
    tentativeTime_ = 0.0;
    currentTimeStepSize_ = 1.0;
    isLocallyConverged_ = false;

    xOld_[0] = p_CC0_;
    x_[0] = p_CC0_;
    p_CA_ = p_CA0_;
    p_CB_ = p_CB0_;
  }

  // non-member ctor
  Teuchos::RCP<pike_test::RxnSingleEq3> 
  rxnSingleEq3(const Teuchos::RCP<pike::MultiphysicsDistributor>& mpd)
  { return Teuchos::rcp(new pike_test::RxnSingleEq3(mpd)); }
  
}
