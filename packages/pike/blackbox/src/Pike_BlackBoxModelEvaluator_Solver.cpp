#include "Pike_BlackBoxModelEvaluator_Solver.hpp"
#include "Pike_Solver.hpp"
#include "Pike_Response.hpp"

namespace pike {

  SolverModelEvaluator::SolverModelEvaluator(const std::string& name) : name_(name) {}

  void SolverModelEvaluator::setSolver(const Teuchos::RCP<pike::Solver>& solver)
  { solver_ = solver; }
  
  Teuchos::RCP<const pike::Solver> SolverModelEvaluator::getSolver() const
  { return solver_; }
  
  Teuchos::RCP<pike::Solver> SolverModelEvaluator::getNonconstSolver() const
  { return solver_; }
  
  std::string SolverModelEvaluator::name() const
  {return name_; }

  bool SolverModelEvaluator::solve()
  { return (solver_->solve() == pike::CONVERGED); }

  bool SolverModelEvaluator::isConverged() const
  {
    return (solver_->getStatus() == pike::CONVERGED);
  }

  bool SolverModelEvaluator::isGloballyConverged() const
  { return true; }

  Teuchos::RCP<const pike::any> SolverModelEvaluator::getResponse(const int i) const
  {
    TEUCHOS_ASSERT(true);
    return Teuchos::null;
  }

  int SolverModelEvaluator::getResponseIndex(const std::string& name) const
  {
    TEUCHOS_ASSERT(true);
    return 0;
  }

  bool SolverModelEvaluator::supportsResponse(const std::string& name) const
  {
    TEUCHOS_ASSERT(true);
    return false;
  }
  
  int SolverModelEvaluator::getNumberOfResponses() const
  {
    TEUCHOS_ASSERT(true);
    return 0;
  }

}
