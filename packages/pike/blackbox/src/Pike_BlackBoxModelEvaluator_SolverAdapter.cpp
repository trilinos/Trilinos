#include "Pike_BlackBoxModelEvaluator_SolverAdapter.hpp"
#include "Pike_Solver.hpp"

namespace pike {

  SolverAdapterModelEvaluator::SolverAdapterModelEvaluator(const std::string& myName) : name_(myName) {}

  void SolverAdapterModelEvaluator::setSolver(const Teuchos::RCP<pike::Solver>& solver)
  {
    solver_ = solver;

    std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > models = 
      solver_->getModelEvaluators();

    parameterNames_.clear();
    parameterNameToIndex_.clear();
    parameterIndexToModelIndices_.clear();
    for (std::size_t m = 0; m < models.size(); ++m) {
      for (int p=0; p < models[m]->getNumberOfParameters(); ++p) {

	// Make sure that different model evaluators don't have the
	// same parameter, otherwise the name to index map will point
	// to only the last parameter. We could change this to support
	// same parameter names in the future, but this is dangerous
	// and possibly not what the user intended.
	std::map<std::string,int>::const_iterator check = parameterNameToIndex_.find(models[m]->getParameterName(p));
	TEUCHOS_TEST_FOR_EXCEPTION(check != parameterNameToIndex_.end(), std::runtime_error,
				   "Error: pike::SolverAdapterModelEvaluator::setSolver() - In the solver adapter \"" 
				   << this->name()
				   << "\", the parameter \"" 
				   << models[m]->getParameterName(p) 
				   << "\" already exists in another model evaluator.  You currently can not have the same parameter name in multiple model evaluators that are both registered to the same SolverAdapterModelEvaluator object!");
	
	parameterNameToIndex_[models[m]->getParameterName(p)] = parameterNames_.size();
	parameterNames_.push_back(models[m]->getParameterName(p));
	parameterIndexToModelIndices_.push_back(std::make_pair(m,p));
      }
    }

    responseNames_.clear();
    responseNameToIndex_.clear();
    responseIndexToModelIndices_.clear();
    for (std::size_t m = 0; m < models.size(); ++m) {
      for (int r=0; r < models[m]->getNumberOfResponses(); ++r) {

	// Make sure that multiple model evaluators do not have the same response
	std::map<std::string,int>::const_iterator check = responseNameToIndex_.find(models[m]->getResponseName(r));
	TEUCHOS_TEST_FOR_EXCEPTION(check != responseNameToIndex_.end(), std::runtime_error,
				   "Error: pike::SolverAdapterModelEvaluator::setSolver() - In the solver adapter \"" 
				   << this->name()
				   << "\", the response \"" 
				   << models[m]->getResponseName(r) 
				   << "\" already exists in another model evaluator.  You can not have the same response name in multiple model evaluators that are both registered to the same SolverAdapterModelEvaluator object!");
	
	responseNameToIndex_[models[m]->getResponseName(r)] = responseNames_.size();
	responseNames_.push_back(models[m]->getResponseName(r));
	responseIndexToModelIndices_.push_back(std::make_pair(m,r));
      }
    }
  }
  
  Teuchos::RCP<const pike::Solver> SolverAdapterModelEvaluator::getSolver() const
  { return solver_; }
  
  Teuchos::RCP<pike::Solver> SolverAdapterModelEvaluator::getNonconstSolver() const
  { return solver_; }
  
  std::string SolverAdapterModelEvaluator::name() const
  {return name_; }

  void SolverAdapterModelEvaluator::solve()
  {
    solver_->reset();
    solver_->solve();
  }

  bool SolverAdapterModelEvaluator::isLocallyConverged() const
  {
    return (solver_->getStatus() == pike::CONVERGED);
  }

  bool SolverAdapterModelEvaluator::isGloballyConverged() const
  { return true; }

  bool SolverAdapterModelEvaluator::supportsParameter(const std::string& pName) const
  {
    return (parameterNameToIndex_.find(pName) !=  parameterNameToIndex_.end());    
  }

  int SolverAdapterModelEvaluator::getNumberOfParameters() const
  {
    return parameterNames_.size();
  }

  std::string SolverAdapterModelEvaluator::getParameterName(const int l) const
  {
    return parameterNames_[l];
  }
  
  int SolverAdapterModelEvaluator::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(parameterNameToIndex_.find(pName) == parameterNameToIndex_.end(),
			       std::logic_error,
			       "The parameter name \"" << pName << "\" does not exist!");
    return parameterNameToIndex_.find(pName)->second;
  }

  void SolverAdapterModelEvaluator::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    // Not ideal.  const_cast or friend class with nonconst private
    // accessor or put public nonconst accessor on solver base.  None
    // are appealing.  This best protects users.
    const_cast<pike::BlackBoxModelEvaluator&>(*(solver_->getModelEvaluators()[responseIndexToModelIndices_[l].first])).setParameter(responseIndexToModelIndices_[l].second,p);
  }

  bool SolverAdapterModelEvaluator::supportsResponse(const std::string& rName) const
  {
    return (responseNameToIndex_.find(rName) !=  responseNameToIndex_.end());
  }
  
  int SolverAdapterModelEvaluator::getNumberOfResponses() const
  {
    return responseNames_.size();
  }

  std::string SolverAdapterModelEvaluator::getResponseName(const int i) const
  {
    return responseNames_[i];
  }

  int SolverAdapterModelEvaluator::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(responseNameToIndex_.find(rName) == responseNameToIndex_.end(),
			       std::logic_error,
			       "The response name \"" << rName << "\" does not exist!");
    return responseNameToIndex_.find(rName)->second;
  }

  Teuchos::ArrayView<const double> SolverAdapterModelEvaluator::getResponse(const int i) const
  {
    return solver_->getModelEvaluators()[responseIndexToModelIndices_[i].first]->getResponse(responseIndexToModelIndices_[i].second);
  }

}
