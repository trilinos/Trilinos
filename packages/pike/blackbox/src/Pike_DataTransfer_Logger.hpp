#ifndef PIKE_DATA_TRANSFER_LOGGER_HPP
#define PIKE_DATA_TRANSFER_LOGGER_HPP

#include "Pike_DataTransfer.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>
#include <string>

namespace pike {
  
  /** \brief A DataTransfer decorator that logs certain method calls.

      Currently, this only logs the doTransfer() method.
   */
  class DataTransferLogger : public pike::DataTransfer {
    
  public:
    
    DataTransferLogger(const Teuchos::RCP<pike::DataTransfer>& transfer);

    void setLog(const Teuchos::RCP<std::vector<std::string> >& log);

    Teuchos::RCP<const std::vector<std::string> > getLog() const;

    Teuchos::RCP<std::vector<std::string> > getNonConstLog() const;

    std::string name() const;

    bool doTransfer(const pike::Solver& solver);

    bool transferSucceeded() const;

    const std::vector<std::string>& getSourceModelNames() const;
    
    const std::vector<std::string>& getTargetModelNames() const;

  private:
    
    Teuchos::RCP<std::vector<std::string> > log_;
    Teuchos::RCP<pike::DataTransfer> transfer_;
  };

  /** Non-member ctor
      \relates DataTransferLogger
  */
  Teuchos::RCP<pike::DataTransferLogger> 
  dataTransferLogger(const Teuchos::RCP<pike::DataTransfer>& transfer);
}

#endif
