#ifndef PIKE_MOCK_DATA_TRANSFER_HPP
#define PIKE_MOCK_DATA_TRANSFER_HPP

#include "Pike_DataTransfer.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos { template<typename T> class Comm; }

namespace pike_test {

  /** \brief Mock data transfer for unit testing

      Users can set the solver iteration to either converge or fail on
      (and can choose whether it is a local or global convergence
      failure).
   */
  class MockDataTransfer : public pike::DataTransfer {

  public:

    MockDataTransfer(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		     const std::string& myName,
		     const std::vector<std::string>& sourceModelNames,
		     const std::vector<std::string>& targetModelNames);
    
    //@{ DataTransfer derived methods
    
    virtual std::string name() const;

    virtual bool doTransfer(const pike::Solver& solver);

    virtual bool transferSucceeded() const;

    virtual const std::vector<std::string>& getSourceModelNames() const;
    
    virtual const std::vector<std::string>& getTargetModelNames() const;
    //@}

  private:
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    std::string name_;
    std::vector<std::string> sourceModelNames_;
    std::vector<std::string> targetModelNames_;
  };

  /** \brief non-member ctor
      \relates MockDataTransfer
  */
  Teuchos::RCP<pike_test::MockDataTransfer> 
  mockDataTransfer(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
		   const std::string& name,
		   const std::vector<std::string>& sourceModelNames,
		   const std::vector<std::string>& targetModelNames);
  

}

#endif
