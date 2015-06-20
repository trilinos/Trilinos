#ifndef PIKE_TEST_MOCK_USER_STATUS_TEST_FACTORY_HPP
#define PIKE_TEST_MOCK_USER_STATUS_TEST_FACTORY_HPP

#include "Pike_StatusTest_AbstractFactory.hpp"

namespace pike_test {

  class UserStatusTestFactory : public pike::StatusTestAbstractFactory {
    
  public:
    UserStatusTestFactory(const std::string& myTestType);

    bool supportsType(const std::string& type) const;
    
    Teuchos::RCP<pike::StatusTest> 
    buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p) const;

  private:
    std::string myType_;
  };

}

#endif
