#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_LOGGER_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_LOGGER_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>
#include <string>

namespace pike {

  /** \brief A BlackBoxModelEvaluator decorator that logs certain method calls.

      Currently, this only logs the solve() and getResponse() methods.
   */
  class ModelEvaluatorLogger : public pike::BlackBoxModelEvaluator {

  public:

    ModelEvaluatorLogger(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model);

    void setLog(const Teuchos::RCP<std::vector<std::string> >& log);

    Teuchos::RCP<const std::vector<std::string> > getLog() const;

    Teuchos::RCP<std::vector<std::string> > getNonConstLog() const;

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
    
    Teuchos::RCP<std::vector<std::string> > log_;
    Teuchos::RCP<pike::BlackBoxModelEvaluator> model_;
  };

  /** \brief Non-member ctor
      \relates ModelEvaluatorLogger
  */
  Teuchos::RCP<ModelEvaluatorLogger>
  modelEvaluatorLogger(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model);

}

#endif
