#include "Pike_Mock_UserStatusTestFactory.hpp"
#include "Pike_StatusTest_MaxIterations.hpp"
namespace pike_test {

  UserStatusTestFactory::UserStatusTestFactory(const std::string& myType)
  {
    myType_ = myType;
  }
  
  bool UserStatusTestFactory::supportsType(const std::string& type) const
  {
    return (type == myType_);
  }
  
  Teuchos::RCP<pike::StatusTest> 
  UserStatusTestFactory::buildStatusTests(const Teuchos::RCP<Teuchos::ParameterList>& p) const
  {
    // Instead of wasting time writing a new user status test, we will
    // reuse a base one here to test the factory.  To reuse, need to
    // change the type so that plist validation works.
    p->set("Type","Maximum Iterations");
    Teuchos::RCP<pike::MaxIterations> test = Teuchos::rcp(new pike::MaxIterations);
    test->setParameterList(p);
    return test;
  }
  
}
