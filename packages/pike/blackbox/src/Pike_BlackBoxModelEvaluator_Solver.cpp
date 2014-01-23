#include "Pike_BlackBoxModelEvaluator_Solver.hpp"
#include "Pike_Solver.hpp"

namespace pike {

  SolverModelEvaluator::SolverModelEvaluator(const std::string& name) : name_(name) {}

  void SolverModelEvaluator::setSolver(const Teuchos::RCP<pike::Solver>& solver)
  {
    solver_ = solver;

    std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> > models = 
      solver_->getModelEvaluators();

    responseNames_.clear();
    responseNameToIndex_.clear();
    responseIndexToModelIndices_.clear();
    typedef std::vector<Teuchos::RCP<const pike::BlackBoxModelEvaluator> >::const_iterator it;
    for (int m = 0; m < models.size(); ++m) {
      for (int r=0; r < models[m]->getNumberOfResponses(); ++r) {
	responseNameToIndex_[models[m]->getResponseName(r)] = responseNames_.size();
	responseNames_.push_back(models[m]->getResponseName(r));
	responseIndexToModelIndices_.push_back(std::make_pair(m,r));
      }
    }
  }
  
  Teuchos::RCP<const pike::Solver> SolverModelEvaluator::getSolver() const
  { return solver_; }
  
  Teuchos::RCP<pike::Solver> SolverModelEvaluator::getNonconstSolver() const
  { return solver_; }
  
  std::string SolverModelEvaluator::name() const
  {return name_; }

  bool SolverModelEvaluator::solve()
  { return (solver_->solve() == pike::CONVERGED); }

  bool SolverModelEvaluator::isLocallyConverged() const
  {
    return (solver_->getStatus() == pike::CONVERGED);
  }

  bool SolverModelEvaluator::isGloballyConverged() const
  { return true; }

  Teuchos::RCP<const pike::any> SolverModelEvaluator::getResponse(const int i) const
  {
    return solver_->getModelEvaluators()[responseIndexToModelIndices_[i].first]->getResponse(responseIndexToModelIndices_[i].second);
  }

  int SolverModelEvaluator::getResponseIndex(const std::string& name) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(responseNameToIndex_.find(name) == responseNameToIndex_.end(),
			       std::logic_error,
			       "The response name \"" << name << "\" does not exist!");
    return responseNameToIndex_.find(name)->second;
  }

  std::string SolverModelEvaluator::getResponseName(const int i) const
  {
    return responseNames_[i];
  }

  bool SolverModelEvaluator::supportsResponse(const std::string& name) const
  {
    return (responseNameToIndex_.find(name) !=  responseNameToIndex_.end());
  }
  
  int SolverModelEvaluator::getNumberOfResponses() const
  {
    return responseNames_.size();
  }

}
