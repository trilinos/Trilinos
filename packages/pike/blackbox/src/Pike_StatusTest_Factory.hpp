#ifndef PIKE_STATUS_TEST_FACTORY
#define PIKE_STATUS_TEST_FACTORY

#include "Teuchos_RCP.hpp"
#include "Pike_StatusTest_AbstractFactory.hpp"
#include <vector>

namespace pike {

  class StatusTestFactory : public pike::StatusTestAbstractFactory {
    
  public:

    StatusTestFactory();

    bool supportsType(const std::string& type) const;
    
    Teuchos::RCP<pike::StatusTest> 
    buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p) const;

    void addFactory(const Teuchos::RCP<pike::StatusTestAbstractFactory>& f);

  private:
    
    //! Contains user added factories for additional types.
    std::vector<Teuchos::RCP<pike::StatusTestAbstractFactory> > userFactories_;

    //! Contains the vaild types this factory supports.
    std::vector<std::string> myTypes_;
  };

}

#endif
