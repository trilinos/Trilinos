#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_LOGGER_DECORATOR_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_LOGGER_DECORATOR_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"
#include "Pike_Response.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>
#include <string>

namespace pike {

  /** \brief A BlackBoxModelEvaluator decorator that logs certain method calls.

      Currently, this only logs the solve() and getResponse() methods.
   */
  class BBMELoggerDecorator : public pike::BlackBoxModelEvaluator {

  public:

    BBMELoggerDecorator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model);

    void setLog(const Teuchos::RCP<std::vector<std::string> >& log);

    Teuchos::RCP<const std::vector<std::string> > getLog() const;

    Teuchos::RCP<std::vector<std::string> > getNonConstLog() const;

    std::string name() const;

    bool solve();

    bool isConverged() const;

    bool isGloballyConverged() const;

    Teuchos::RCP<pike::Response> getResponse(const int i) const;

    int getResponseIndex(const std::string name) const;

    bool supportsResponse(const std::string name) const;

  private:
    
    Teuchos::RCP<std::vector<std::string> > log_;
    Teuchos::RCP<pike::BlackBoxModelEvaluator> model_;
  };

  /** \brief Non-member ctor
      \relates BBMELoggerDecorator
  */
  Teuchos::RCP<BBMELoggerDecorator>
  bbmeLoggerDecorator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model);

}

#endif
