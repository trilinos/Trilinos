#ifndef PIKE_LINEAR_HEAT_CONDUTION_MODEL_EVALUATOR_HPP
#define PIKE_LINEAR_HEAT_CONDUTION_MODEL_EVALUATOR_HPP

#include "Pike_LinearHeatConductionModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"

namespace pike_test {
  
  LinearHeatConductionModelEvaluator::LinearHeatConductionModelEvaluator(Teuchos::RCP<Teuchos::Comm<int> > comm,
									 std::string name,
									 Mode mode)
    : comm_(comm),
      name_(name),
      mode_(mode),
      k_(1.0),
      q_(1.0),
      T_left_(1.0),
      T_right_(1.0),
      responseValue_(1) // only one response
  {
    if (mode_ == T_RIGHT_IS_RESPONSE) {
      responseMap_["T_right"] = 0;
    }
    else if (mode_ == Q_IS_RESPONSE) {
      responseMap_["q"] = 0;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error the mode is not valid!");
    }
  }
  
  LinearHeatConductionModelEvaluator::~LinearHeatConductionModelEvaluator() {}
  
  const std::string LinearHeatConductionModelEvaluator::name() const
  { return name_; }
  
  void LinearHeatConductionModelEvaluator::solve()
  {
    // solution: T_left - T_right = q / k
    
    if (mode_ == T_RIGHT_IS_RESPONSE) {
      responseValue_[0]->set(T_left_ - q_ / k_);
    }
    else if (mode_ == Q_IS_RESPONSE) {
      responseValue_[0]->set((T_left_ - T_right_) * k_);
    }
  }
  
  bool LinearHeatConductionModelEvaluator::isConverged() const
  { return true; }
  
  Teuchos::RCP<pike::Response> LinearHeatConductionModelEvaluator::getResponse(const int i) const
  {
    return responseValue_[i];
  }
  
  int LinearHeatConductionModelEvaluator::getResponseIndex(const std::string name) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(responseMap_.find(name) == responseMap_.end(),
			       std::logic_error,
			       "Response name is not valid!");
    return responseMap_.find(name)->second;
  }
  
  bool LinearHeatConductionModelEvaluator::supportsResponse(const std::string name) const
  {
    return (responseMap_.find(name) != responseMap_.end());
  }
  
  void LinearHeatConductionModelEvaluator::set_k(const double& k)
  { k_ = k; }

  void LinearHeatConductionModelEvaluator::set_q(const double& q)
  {
    TEUCHOS_ASSERT(mode_ == T_RIGHT_IS_RESPONSE);
    q_ = q;
  }
  
  void LinearHeatConductionModelEvaluator::set_T_left(const double& T_left)
  { T_left_ = T_left; }
  
  void LinearHeatConductionModelEvaluator::set_T_right(const double& T_right)
  {
    TEUCHOS_ASSERT(mode_ == Q_IS_RESPONSE);
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
  linearHeatConductionModelEvaluator(Teuchos::RCP<Teuchos::Comm<int> > comm,
				     std::string name,
				     pike_test::LinearHeatConductionModelEvaluator::Mode mode)
  {
    return Teuchos::rcp(new pike_test::LinearHeatConductionModelEvaluator(comm,name,mode));
  }

}

#endif
