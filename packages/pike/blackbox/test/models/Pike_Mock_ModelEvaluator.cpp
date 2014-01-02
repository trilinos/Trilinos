#ifndef PIKE_MOCK_MODEL_EVALUATOR_HPP
#define PIKE_MOCK_MODEL_EVALUATOR_HPP

#include "Pike_Mock_ModelEvaluator.hpp"
#include "Pike_Solver.hpp"
#include "Teuchos_Assert.hpp"

namespace pike_test {
  
  MockModelEvaluator::MockModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
					 const std::string& name,
					 const Mode mode,
					 const int iterationTrigger,
					 const int responseFreezeIteration)
    : comm_(comm),
      name_(name),
      mode_(mode),
      iterationTrigger_(iterationTrigger),
      responseFreezeIteration_(responseFreezeIteration),
      responseValue_(1) // only one response
  {
    responseMap_["Mock Response"] = 0;
    responseValue_[0] = Teuchos::rcp(new pike::ScalarResponse<double>("Mock Response"));
  }

  std::string MockModelEvaluator::name() const
  { return name_; }
  
  bool MockModelEvaluator::solve()
  {
    TEUCHOS_ASSERT(nonnull(solver_));

    if (responseFreezeIteration_ <= (solver_->getNumberOfIterations()+1) ) {
      // do nothing - freeze the response
    }
    else
      responseValue_[0]->set(responseValue_[0]->get() + 1.0);

    if (iterationTrigger_ == (solver_->getNumberOfIterations()+1) )
      if (mode_ == LOCAL_FAILURE)
	return false;
    
    return true;
  }
  
  bool MockModelEvaluator::isConverged() const
  { 
    TEUCHOS_ASSERT(nonnull(solver_));
    
    if (iterationTrigger_ == solver_->getNumberOfIterations() )
      if (mode_ == LOCAL_FAILURE)
	return false;
  
    return true;
  }

  bool MockModelEvaluator::isGloballyConverged() const
  {
    TEUCHOS_ASSERT(nonnull(solver_));
    
    if (iterationTrigger_ == solver_->getNumberOfIterations() )
      if (mode_ == GLOBAL_CONVERGENCE)
	return true;

    return false;
  }
  
  Teuchos::RCP<pike::Response> MockModelEvaluator::getResponse(const int i) const
  {
    return responseValue_[i];
  }
  
  int MockModelEvaluator::getResponseIndex(const std::string name) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(responseMap_.find(name) == responseMap_.end(),
			       std::logic_error,
			       "Response name is not valid!");
    return responseMap_.find(name)->second;
  }
  
  bool MockModelEvaluator::supportsResponse(const std::string name) const
  {
    return (responseMap_.find(name) != responseMap_.end());
  }

  void MockModelEvaluator::setSolver(const Teuchos::RCP<pike::Solver>& solver)
  {
    solver_ = solver;
  }

  // non-member ctor
  Teuchos::RCP<pike_test::MockModelEvaluator> 
  mockModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		     const std::string& name,
		     const pike_test::MockModelEvaluator::Mode mode,
		     const int iterationTrigger,
		     const int responseFreezeIteration)
  {
    return Teuchos::rcp(new pike_test::MockModelEvaluator(comm,name,mode,iterationTrigger,responseFreezeIteration));
  }

}

#endif
