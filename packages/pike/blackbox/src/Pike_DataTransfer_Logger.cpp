#include "Pike_DataTransfer_Logger.hpp"

namespace pike {

  DataTransferLogger::DataTransferLogger(const Teuchos::RCP<pike::DataTransfer>& transfer)
    : transfer_(transfer)
  {
    log_ = Teuchos::rcp(new std::vector<std::string>);
  }
  
  void DataTransferLogger::setLog(const Teuchos::RCP<std::vector<std::string> >& log)
  {
    log_ = log;
  }
  
  Teuchos::RCP<const std::vector<std::string> > DataTransferLogger::getLog() const
  {
    return log_;
  }
  
  Teuchos::RCP<std::vector<std::string> > DataTransferLogger::getNonConstLog() const
  {
    return log_;
  }
  
  std::string DataTransferLogger::name() const
  {
    return transfer_->name();
  }
  
  bool DataTransferLogger::doTransfer(const pike::Solver& solver)
  {
    log_->push_back(this->name()+": doTransfer()");
    return transfer_->doTransfer(solver);
  }
  
  bool DataTransferLogger::transferSucceeded() const
  {
    return transfer_->transferSucceeded();
  }

  const std::vector<std::string>& DataTransferLogger::getSourceModelNames() const
  {
    return transfer_->getSourceModelNames();
  }
  
  const std::vector<std::string>& DataTransferLogger::getTargetModelNames() const
  {
    return transfer_->getTargetModelNames();
  }

  // Non-member ctor
  Teuchos::RCP<pike::DataTransferLogger> 
  dataTransferLogger(const Teuchos::RCP<pike::DataTransfer>& transfer)
  {
    return Teuchos::rcp(new pike::DataTransferLogger(transfer));
  }
  
}
