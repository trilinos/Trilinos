#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

namespace {

  using Teuchos::rcp;
  using Teuchos::RCP;
  using MueLu::Needs;

  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST(Needs, Constructor)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Needs> needs = rcp(new Needs() );
    TEUCHOS_TEST_EQUALITY(needs != Teuchos::null, true, out, success);
  }

  TEUCHOS_UNIT_TEST(Needs, NeedRequested)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    TEUCHOS_TEST_EQUALITY(needs.IsRequested(aNeed), false, out, success);
    needs.Request(aNeed);
    TEUCHOS_TEST_EQUALITY(needs.IsRequested(aNeed), true, out, success);
  }

  TEUCHOS_UNIT_TEST(Needs, ValueSaved)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    TEUCHOS_TEST_EQUALITY(needs.IsSaved(aNeed), false, out, success);
    needs.Save(aNeed,42);
    TEUCHOS_TEST_EQUALITY(needs.IsSaved(aNeed), true, out, success);
  }

  TEUCHOS_UNIT_TEST(Needs, NumRequests_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    TEST_THROW( needs.NumRequests("nonExistentNeed"), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, NumRequests)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    needs.Request(aNeed);
    needs.Request(aNeed);
    TEUCHOS_TEST_EQUALITY(needs.NumRequests(aNeed), 2, out, success);
  }

  TEUCHOS_UNIT_TEST(Needs, Examine_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    double value=0;
    TEST_THROW( needs.Examine("nonExistentNeed",value), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, SaveAndExamine)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Save(aNeed,trueValue);
    double expectedValue = 0;
    needs.Examine(aNeed,expectedValue);
    TEUCHOS_TEST_EQUALITY(trueValue,expectedValue, out, success);
  }

  TEUCHOS_UNIT_TEST(Needs, CheckOut_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    double value;
    TEST_THROW( needs.CheckOut("nonExistentNeed",value), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, Checkout_Without_Request)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Save(aNeed,trueValue);
    double expectedValue = 0;
    TEST_THROW( needs.CheckOut(aNeed,expectedValue), MueLu::Exceptions::RuntimeError );
  }

  TEUCHOS_UNIT_TEST(Needs, CheckOut)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Save(aNeed,trueValue);
    needs.Request(aNeed);
    needs.Request(aNeed);
    double value = 0;
    needs.CheckOut(aNeed,value);
    TEUCHOS_TEST_EQUALITY(trueValue,value, out, success);
    TEUCHOS_TEST_EQUALITY(needs.NumRequests(aNeed),1, out, success);
    value = 0;
    needs.CheckOut(aNeed,value);
    //try to get the need one too many times
    TEST_THROW( needs.Examine(aNeed,value), std::logic_error );
  }

  class foobarClass {
  public:
    foobarClass() {}
    virtual ~foobarClass() {}
  };

  TEUCHOS_UNIT_TEST(Needs, nonPOD)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    RCP<foobarClass> trueValue = rcp(new foobarClass);
    needs.Save("foobar",trueValue);
    RCP<foobarClass> value;
    needs.Examine("foobar",value);
    TEUCHOS_TEST_EQUALITY(trueValue,value, out, success);
  }

}//namespace <anonymous>
