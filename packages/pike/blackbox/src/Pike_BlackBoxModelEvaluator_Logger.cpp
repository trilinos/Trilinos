
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

  void ModelEvaluatorLogger::solve()
  {
    log_->push_back(this->name()+": solve()");
    model_->solve();
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
  
  int ModelEvaluatorLogger::getResponseIndex(const std::string& rName) const
  {
    return model_->getResponseIndex(rName);
  }

  std::string ModelEvaluatorLogger::getResponseName(const int i) const
  {
    return model_->getResponseName(i);
  }

  bool ModelEvaluatorLogger::supportsResponse(const std::string& rName) const
  {
    return model_->supportsResponse(rName);
  }

  int ModelEvaluatorLogger::getNumberOfResponses() const
  {
    return model_->getNumberOfResponses();
  }

  bool ModelEvaluatorLogger::supportsParameter(const std::string& pName) const
  {
    return model_->supportsParameter(pName);
  }

  int ModelEvaluatorLogger::getNumberOfParameters() const
  {
    return model_->getNumberOfParameters();
  }

  std::string ModelEvaluatorLogger::getParameterName(const int l) const
  {
    return model_->getParameterName(l);
  }

  int ModelEvaluatorLogger::getParameterIndex(const std::string& pName) const
  {
    return model_->getParameterIndex(pName);
  }

  void ModelEvaluatorLogger::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    log_->push_back(this->name()+": setParameter(l,p)");
    return model_->setParameter(l,p);
  }

  bool ModelEvaluatorLogger::isTransient() const
  {
    return model_->isTransient();
  }

  double ModelEvaluatorLogger::getCurrentTime() const
  {
    return model_->getCurrentTime();
  }

  double ModelEvaluatorLogger::getTentativeTime() const
  {
    return model_->getTentativeTime();
  }

  bool ModelEvaluatorLogger::solvedTentativeStep() const
  {
    return model_->solvedTentativeStep();
  }

  double ModelEvaluatorLogger::getCurrentTimeStepSize() const
  {
    return model_->getCurrentTimeStepSize();
  }

  double ModelEvaluatorLogger::getDesiredTimeStepSize() const
  {
    return model_->getDesiredTimeStepSize();
  }

  double ModelEvaluatorLogger::getMaxTimeStepSize() const
  {
    return model_->getMaxTimeStepSize();
  }

  void ModelEvaluatorLogger::setNextTimeStepSize(const double& dt)
  {
    return model_->setNextTimeStepSize(dt);
  }

  void ModelEvaluatorLogger::acceptTimeStep()
  {
    log_->push_back(this->name()+": acceptTimeStep()");
    return model_->acceptTimeStep();
  }

  //! Non-member ctor
  Teuchos::RCP<ModelEvaluatorLogger>
  modelEvaluatorLogger(const Teuchos::RCP<pike::BlackBoxModelEvaluator>& model)
  {
    return Teuchos::rcp(new ModelEvaluatorLogger(model));
  }

}
