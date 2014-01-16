
#include "Pike_BlackBoxModelEvaluator_Logger.hpp"

namespace pike {

  ModelEvaluatorLogger::ModelEvaluatorLogger(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
    : model_(model)
  { 
    log_ = Teuchos::rcp(new std::vector<std::string>);
  }
  
  void ModelEvaluatorLogger::setLog(const Teuchos::RCP<std::vector<std::string> >& log)
  {
    log_ = log;
  }

  Teuchos::RCP<const std::vector<std::string> > ModelEvaluatorLogger::getLog() const
  {
    return log_;
  }

  Teuchos::RCP<std::vector<std::string> > ModelEvaluatorLogger::getNonConstLog() const
  {
    return log_;
  }

  std::string ModelEvaluatorLogger::name() const
  {
    return model_->name();
  }

  bool ModelEvaluatorLogger::solve()
  {
    log_->push_back(this->name()+": solve()");
    return model_->solve();
  }

  bool ModelEvaluatorLogger::isConverged() const
  { 
    return model_->isConverged();
  }
  
  bool ModelEvaluatorLogger::isGloballyConverged() const
  {
    return model_->isGloballyConverged();
  }
  
  Teuchos::RCP<const pike::any> ModelEvaluatorLogger::getResponse(const int i) const
  {
    log_->push_back(this->name()+": getResponse()");
    return model_->getResponse(i);
  }
  
  int ModelEvaluatorLogger::getResponseIndex(const std::string& name) const
  {
    return model_->getResponseIndex(name);
  }

  std::string ModelEvaluatorLogger::getResponseName(const int i) const
  {
    return model_->getResponseName(i);
  }

  bool ModelEvaluatorLogger::supportsResponse(const std::string& name) const
  {
    return model_->supportsResponse(name);
  }

  int ModelEvaluatorLogger::getNumberOfResponses() const
  {
    return model_->getNumberOfResponses();
  }

  //! Non-member ctor
  Teuchos::RCP<ModelEvaluatorLogger>
  modelEvaluatorLogger(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
  {
    return Teuchos::rcp(new ModelEvaluatorLogger(model));
  }

}
