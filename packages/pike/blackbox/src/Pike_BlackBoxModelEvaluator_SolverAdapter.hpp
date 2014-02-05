#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_SOLVER_ADAPTER_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_SOLVER_ADAPTER_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_RCP.hpp"
#include <utility>
#include <string>
#include <map>

namespace pike {

  class Solver;
  class Response;

  //! Decorator to represent a pike::Solver as a BlackBoxModelEvaluator for hierarchical solves.
  class SolverAdapterModelEvaluator : public pike::BlackBoxModelEvaluator {

  public:

    SolverAdapterModelEvaluator(const std::string& myName);

    void setSolver(const Teuchos::RCP<pike::Solver>& solver);

    Teuchos::RCP<const pike::Solver> getSolver() const;

    Teuchos::RCP<pike::Solver> getNonconstSolver() const;

    // Derived from base
    std::string name() const;
    void solve();
    bool isLocallyConverged() const;
    bool isGloballyConverged() const;
    bool supportsParameter(const std::string& pName) const;
    int getNumberOfParameters() const;
    std::string getParameterName(const int l) const;
    int getParameterIndex(const std::string& pName) const;
    void setParameter(const int l, const Teuchos::ArrayView<const double>& p);
    bool supportsResponse(const std::string& rName) const;
    int getNumberOfResponses() const;
    std::string getResponseName(const int i) const;
    int getResponseIndex(const std::string& rName) const;
    Teuchos::ArrayView<const double> getResponse(const int i) const;

  private:
    std::string name_;
    Teuchos::RCP<pike::Solver> solver_;

    std::map<std::string,int> parameterNameToIndex_;
    std::vector<std::string> parameterNames_;
    //! Stores the model index and the response index in that model for the response. 
    std::vector<std::pair<int,int> > parameterIndexToModelIndices_;

    std::map<std::string,int> responseNameToIndex_;
    std::vector<std::string> responseNames_;
    //! Stores the model index and the response index in that model for the response. 
    std::vector<std::pair<int,int> > responseIndexToModelIndices_;
  };

}

#endif
