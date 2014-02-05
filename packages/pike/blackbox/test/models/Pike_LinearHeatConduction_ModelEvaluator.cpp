#include "Pike_LinearHeatConduction_ModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

namespace pike_test {
  
  LinearHeatConductionModelEvaluator::LinearHeatConductionModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
									 const std::string& myName,
									 const Mode mode)
    : comm_(comm),
      name_(myName),
      mode_(mode),
      k_(1.0),
      q_(1.0),
      T_left_(1.0),
      T_right_(1.0),
      responseNames_(1), // only one response
      responseValues_(1)
  {

    if (mode_ == T_RIGHT_IS_RESPONSE) {
      responseMap_["T_right"] = 0;
      responseNames_[0] = "T_right";
      responseValues_[0].resize(1);
      responseValues_[0][0] = T_right_;
    }
    else if (mode_ == Q_IS_RESPONSE) {
      responseMap_["q"] = 0;
      responseNames_[0] = "q";
      responseValues_[0].resize(1);
      responseValues_[0][0] = q_;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error the mode is not valid!");
    }
  }

  std::string LinearHeatConductionModelEvaluator::name() const
  { return name_; }
  
  void LinearHeatConductionModelEvaluator::solve()
  {
    // solution: T_left - T_right = q / k
    
    if (mode_ == T_RIGHT_IS_RESPONSE) {
      T_right_ = T_left_ - q_ / k_;
      responseValues_[0][0] = T_right_;
    }
    else if (mode_ == Q_IS_RESPONSE) {
      q_ = (T_left_ - T_right_) * k_;
      responseValues_[0][0] = q_;
    }
  }
  
  bool LinearHeatConductionModelEvaluator::isLocallyConverged() const
  { return true; }

  bool LinearHeatConductionModelEvaluator::isGloballyConverged() const
  { return true; }
  
  Teuchos::ArrayView<const double> LinearHeatConductionModelEvaluator::getResponse(const int i) const
  {
    return Teuchos::ArrayView<const double>(responseValues_[i]);
  }
  
  int LinearHeatConductionModelEvaluator::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(responseMap_.find(rName) == responseMap_.end(),
			       std::logic_error,
			       "Response name \"" << rName << "\"is not valid!");
    return responseMap_.find(rName)->second;
  }
  
  std::string LinearHeatConductionModelEvaluator::getResponseName(const int i) const
  {
    TEUCHOS_ASSERT( (i>=0) && (i<Teuchos::as<int>(responseNames_.size())) );
    return responseNames_[i];
  }

  bool LinearHeatConductionModelEvaluator::supportsResponse(const std::string& rName) const
  {
    return (responseMap_.find(rName) != responseMap_.end());
  }
  
  int LinearHeatConductionModelEvaluator::getNumberOfResponses() const
  {
    return Teuchos::as<int>(responseMap_.size());
  }

  void LinearHeatConductionModelEvaluator::set_k(const double& k)
  { k_ = k; }

  void LinearHeatConductionModelEvaluator::set_q(const double& q)
  {
    q_ = q;
  }
  
  void LinearHeatConductionModelEvaluator::set_T_left(const double& T_left)
  { T_left_ = T_left; }
  
  void LinearHeatConductionModelEvaluator::set_T_right(const double& T_right)
  {
    T_right_ = T_right;
  }
  
  double LinearHeatConductionModelEvaluator::get_q() const
  { return q_; }
  
  double LinearHeatConductionModelEvaluator::get_K() const
  { return k_; }

  double LinearHeatConductionModelEvaluator::get_T_left() const
  { return T_left_; }
  
  double LinearHeatConductionModelEvaluator::get_T_right() const
  { return T_right_; }
  
  // non-member ctor
  Teuchos::RCP<pike_test::LinearHeatConductionModelEvaluator> 
  linearHeatConductionModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
				     const std::string& name,
				     const pike_test::LinearHeatConductionModelEvaluator::Mode mode)
  {
    return Teuchos::rcp(new pike_test::LinearHeatConductionModelEvaluator(comm,name,mode));
  }

}
