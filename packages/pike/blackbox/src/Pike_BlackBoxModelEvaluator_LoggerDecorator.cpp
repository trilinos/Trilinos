
#include "Pike_BlackBoxModelEvaluator_LoggerDecorator.hpp"

namespace pike {

  ModelLoggerDecorator::ModelLoggerDecorator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
    : model_(model)
  { 
    log_ = Teuchos::rcp(new std::vector<std::string>);
  }
  
  void ModelLoggerDecorator::setLog(const Teuchos::RCP<std::vector<std::string> >& log)
  {
    log_ = log;
  }

  Teuchos::RCP<const std::vector<std::string> > ModelLoggerDecorator::getLog() const
  {
    return log_;
  }

  Teuchos::RCP<std::vector<std::string> > ModelLoggerDecorator::getNonConstLog() const
  {
    return log_;
  }

  std::string ModelLoggerDecorator::name() const
  {
    return model_->name();
  }

  bool ModelLoggerDecorator::solve()
  {
    log_->push_back(this->name()+": solve()");
    return model_->solve();
  }

  bool ModelLoggerDecorator::isConverged() const
  { 
    return model_->isConverged();
  }
  
  bool ModelLoggerDecorator::isGloballyConverged() const
  {
    return model_->isGloballyConverged();
  }
  
  Teuchos::RCP<pike::Response> ModelLoggerDecorator::getResponse(const int i) const
  {
    log_->push_back(this->name()+": getResponse()");
    return model_->getResponse(i);
  }
  
  int ModelLoggerDecorator::getResponseIndex(const std::string name) const
  {
    return model_->getResponseIndex(name);
  }

  bool ModelLoggerDecorator::supportsResponse(const std::string name) const
  {
    return model_->supportsResponse(name);
  }

  //! Non-member ctor
  Teuchos::RCP<ModelLoggerDecorator>
  modelLoggerDecorator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
  {
    return Teuchos::rcp(new ModelLoggerDecorator(model));
  }

}
