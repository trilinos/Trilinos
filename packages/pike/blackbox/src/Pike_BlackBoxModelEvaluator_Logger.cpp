
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

  bool ModelEvaluatorLogger::isLocallyConverged() const
  { 
    return model_->isLocallyConverged();
  }
  
  bool ModelEvaluatorLogger::isGloballyConverged() const
  {
    return model_->isGloballyConverged();
  }
  
  Teuchos::ArrayView<const double> ModelEvaluatorLogger::getResponse(const int i) const
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

  bool ModelEvaluatorLogger::supportsParameter(const std::string& name) const
  {
    return model_->supportsParameter(name);
  }

  int ModelEvaluatorLogger::getNumberOfParameters() const
  {
    return model_->getNumberOfParameters();
  }

  std::string ModelEvaluatorLogger::getParameterName(const int l) const
  {
    return model_->getParameterName(l);
  }

  int ModelEvaluatorLogger::getParameterIndex(const std::string& name) const
  {
    return model_->getParameterIndex(name);
  }

  void ModelEvaluatorLogger::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    return model_->setParameter(l,p);
  }

  //! Non-member ctor
  Teuchos::RCP<ModelEvaluatorLogger>
  modelEvaluatorLogger(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
  {
    return Teuchos::rcp(new ModelEvaluatorLogger(model));
  }

}
