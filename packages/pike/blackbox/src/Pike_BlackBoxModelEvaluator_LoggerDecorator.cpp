
#include "Pike_BlackBoxModelEvaluator_LoggerDecorator.hpp"

namespace pike {

  BBMELoggerDecorator::BBMELoggerDecorator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
    : model_(model)
  { 
    log_ = Teuchos::rcp(new std::vector<std::string>);
  }
  
  void BBMELoggerDecorator::setLog(const Teuchos::RCP<std::vector<std::string> >& log)
  {
    log_ = log;
  }

  Teuchos::RCP<const std::vector<std::string> > BBMELoggerDecorator::getLog() const
  {
    return log_;
  }

  Teuchos::RCP<std::vector<std::string> > BBMELoggerDecorator::getNonConstLog() const
  {
    return log_;
  }

  std::string BBMELoggerDecorator::name() const
  {
    return model_->name();
  }

  bool BBMELoggerDecorator::solve()
  {
    log_->push_back(this->name()+": solve()");
    return model_->solve();
  }

  bool BBMELoggerDecorator::isConverged() const
  { 
    return model_->isConverged();
  }
  
  bool BBMELoggerDecorator::isGloballyConverged() const
  {
    return model_->isGloballyConverged();
  }
  
  Teuchos::RCP<pike::Response> BBMELoggerDecorator::getResponse(const int i) const
  {
    log_->push_back(this->name()+": getResponse()");
    return model_->getResponse(i);
  }
  
  int BBMELoggerDecorator::getResponseIndex(const std::string name) const
  {
    return model_->getResponseIndex(name);
  }

  bool BBMELoggerDecorator::supportsResponse(const std::string name) const
  {
    return model_->supportsResponse(name);
  }

  //! Non-member ctor
  Teuchos::RCP<BBMELoggerDecorator>
  bbmeLoggerDecorator(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
  {
    return Teuchos::rcp(new BBMELoggerDecorator(model));
  }

}
