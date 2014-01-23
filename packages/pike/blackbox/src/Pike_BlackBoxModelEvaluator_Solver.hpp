#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_SOLVER_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_SOLVER_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_RCP.hpp"
#include <utility>
#include <string>
#include <map>

namespace pike {

  class Solver;
  class Response;

  //! Decorator to represent a pike::Solver as a BlackBoxModelEvaluator for hierarchical solves.
  class SolverModelEvaluator : pike::BlackBoxModelEvaluator {

  public:

    SolverModelEvaluator(const std::string& name);

    void setSolver(const Teuchos::RCP<pike::Solver>& solver);

    Teuchos::RCP<const pike::Solver> getSolver() const;

    Teuchos::RCP<pike::Solver> getNonconstSolver() const;

    // Derived from base
    std::string name() const;
    bool solve();
    bool isLocallyConverged() const;
    bool isGloballyConverged() const;
    Teuchos::RCP<const pike::any> getResponse(const int i) const;
    int getResponseIndex(const std::string& name) const;
    std::string getResponseName(const int i) const;
    bool supportsResponse(const std::string& name) const;
    int getNumberOfResponses() const;

  private:
    std::string name_;
    Teuchos::RCP<pike::Solver> solver_;
    std::map<std::string,int> responseNameToIndex_;
    std::vector<std::string> responseNames_;
    //! Stores the model index and the response index in that model for the response. 
    std::vector<std::pair<int,int> > responseIndexToModelIndices_;
  };

}

#endif
