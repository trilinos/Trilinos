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
      parameterNames_(1),
      responseNames_(1), // only one response
      responseValues_(1)
  {

    if (mode_ == T_RIGHT_IS_RESPONSE) {
      parameterMap_["q"] = 0;
      parameterNames_[0] = "q";
      responseMap_["T_right"] = 0;
      responseNames_[0] = "T_right";
      responseValues_[0].resize(1);
      responseValues_[0][0] = T_right_;
    }
    else if (mode_ == Q_IS_RESPONSE) {
      parameterMap_["T_right"] = 0;
      parameterNames_[0] = "T_right";
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
  
  bool LinearHeatConductionModelEvaluator::supportsParameter(const std::string& pName) const
  {
    return (parameterMap_.find(pName) != parameterMap_.end());
  }

  int LinearHeatConductionModelEvaluator::getNumberOfParameters() const
  {
    return parameterMap_.size();
  }

  std::string LinearHeatConductionModelEvaluator::getParameterName(const int l) const
  {
    TEUCHOS_ASSERT( (l>=0) && (l<Teuchos::as<int>(parameterNames_.size())) );
    return parameterNames_[l];
  }

  int LinearHeatConductionModelEvaluator::getParameterIndex(const std::string& pName) const
  {
    std::map<std::string,int>::const_iterator l = parameterMap_.find(pName);
    TEUCHOS_TEST_FOR_EXCEPTION(l == parameterMap_.end(),
			       std::logic_error,
			       "Parameter name \"" << pName << "\"is not valid!");
    return l->second;
  }

  void LinearHeatConductionModelEvaluator::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_ASSERT( (l>=0) && (l<Teuchos::as<int>(parameterNames_.size())) );
    if (mode_ == T_RIGHT_IS_RESPONSE)
      this->set_q(p[0]);
    else
      this->set_T_right(p[0]);
  }

  Teuchos::ArrayView<const double> LinearHeatConductionModelEvaluator::getResponse(const int i) const
  {
    return Teuchos::ArrayView<const double>(responseValues_[i]);
  }
  
  int LinearHeatConductionModelEvaluator::getResponseIndex(const std::string& rName) const
  {
    std::map<std::string,int>::const_iterator  i = responseMap_.find(rName);
    TEUCHOS_TEST_FOR_EXCEPTION(i == responseMap_.end(),
			       std::logic_error,
			       "Response name \"" << rName << "\"is not valid!");
    return i->second;
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
