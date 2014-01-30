#include "Pike_Mock_ModelEvaluator.hpp"
#include "Pike_Solver.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

namespace pike_test {
  
  MockModelEvaluator::MockModelEvaluator(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
					 const std::string& myName,
					 const Mode mode,
					 const int iterationTrigger,
					 const int responseFreezeIteration)
    : comm_(comm),
      name_(myName),
      mode_(mode),
      iterationTrigger_(iterationTrigger),
      responseFreezeIteration_(responseFreezeIteration),
      responseNames_(1), // only one response
      responseValues_(1)
  {
    responseMap_["Mock Response"] = 0;
    responseNames_[0] = "Mock Response";
    responseValues_[0].resize(1);
    responseValues_[0][0] = 0.0;
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
      responseValues_[0][0] += 1.0;

    if (iterationTrigger_ == (solver_->getNumberOfIterations()+1) )
      if (mode_ == LOCAL_FAILURE)
	return false;
    
    return true;
  }
  
  bool MockModelEvaluator::isLocallyConverged() const
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
  
  Teuchos::ArrayView<const double> MockModelEvaluator::getResponse(const int i) const
  {
    return Teuchos::ArrayView<const double>(responseValues_[i]);
  }
  
  int MockModelEvaluator::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(responseMap_.find(rName) == responseMap_.end(),
			       std::logic_error,
			       "Response name \"" << rName << "\" is not valid!");
    return responseMap_.find(rName)->second;
  }
  
  std::string MockModelEvaluator::getResponseName(const int i) const
  {
    TEUCHOS_ASSERT( (i >=0) && (i<Teuchos::as<int>(responseNames_.size())) );
    return responseNames_[i];
  }

  bool MockModelEvaluator::supportsResponse(const std::string& rName) const
  {
    return (responseMap_.find(rName) != responseMap_.end());
  }

  int MockModelEvaluator::getNumberOfResponses() const
  {
    return Teuchos::as<int>(responseMap_.size());
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
