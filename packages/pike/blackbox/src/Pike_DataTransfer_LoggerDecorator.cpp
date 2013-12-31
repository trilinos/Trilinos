#include "Pike_DataTransfer_LoggerDecorator.hpp"

namespace pike {

  DataTransferLoggerDecorator::DataTransferLoggerDecorator(const Teuchos::RCP<pike::DataTransfer>& transfer)
    : transfer_(transfer)
  {
    log_ = Teuchos::rcp(new std::vector<std::string>);
  }
  
  void DataTransferLoggerDecorator::setLog(const Teuchos::RCP<std::vector<std::string> >& log)
  {
    log_ = log;
  }
  
  Teuchos::RCP<const std::vector<std::string> > DataTransferLoggerDecorator::getLog() const
  {
    return log_;
  }
  
  Teuchos::RCP<std::vector<std::string> > DataTransferLoggerDecorator::getNonConstLog() const
  {
    return log_;
  }
  
  std::string DataTransferLoggerDecorator::name() const
  {
    return transfer_->name();
  }
  
  bool DataTransferLoggerDecorator::doTransfer(const pike::Solver& solver)
  {
    log_->push_back(this->name()+": doTransfer()");
    return transfer_->doTransfer(solver);
  }
  
  bool DataTransferLoggerDecorator::transferSucceeded() const
  {
    return transfer_->transferSucceeded();
  }

  const std::vector<std::string>& DataTransferLoggerDecorator::getSourceModelNames() const
  {
    return transfer_->getSourceModelNames();
  }
  
  const std::vector<std::string>& DataTransferLoggerDecorator::getTargetModelNames() const
  {
    return transfer_->getTargetModelNames();
  }

  // Non-member ctor
  Teuchos::RCP<pike::DataTransferLoggerDecorator> 
  dataTransferLoggerDecorator(const Teuchos::RCP<pike::DataTransfer>& transfer)
  {
    return Teuchos::rcp(new pike::DataTransferLoggerDecorator(transfer));
  }
  
}
