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

  TEUCHOS_UNIT_TEST(Needs, ValueIsAvailable)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    TEUCHOS_TEST_EQUALITY(needs.IsAvailable(aNeed), false, out, success);
    needs.Set(aNeed,42);
    TEUCHOS_TEST_EQUALITY(needs.IsAvailable(aNeed), true, out, success);
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

  TEUCHOS_UNIT_TEST(Needs, Get_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    double value=0;
    TEST_THROW( needs.Get("nonExistentNeed",value), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, SetAndGet)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Set(aNeed,trueValue);
    double expectedValue = 0;
    needs.Get(aNeed,expectedValue);
    TEUCHOS_TEST_EQUALITY(trueValue,expectedValue, out, success);
  }

  TEUCHOS_UNIT_TEST(Needs, Release_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    TEST_THROW( needs.Release("nonExistentNeed"), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, Release_Without_Request)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Set(aNeed,trueValue);
    double expectedValue = 0;
    TEST_THROW( needs.Get(aNeed,expectedValue), MueLu::Exceptions::RuntimeError );
    TEST_THROW( needs.Release(aNeed), MueLu::Exceptions::RuntimeError );
  }

  TEUCHOS_UNIT_TEST(Needs, Release)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Set(aNeed,trueValue);
    needs.Request(aNeed);
    needs.Request(aNeed);
    double value = 0;
    needs.Get(aNeed,value);
    needs.Release(aNeed);
    TEUCHOS_TEST_EQUALITY(trueValue,value, out, success);
    TEUCHOS_TEST_EQUALITY(needs.NumRequests(aNeed),1, out, success);
    value = 0;
    needs.Get(aNeed,value);
    needs.Release(aNeed);
    //try to get the need one too many times
    //JG TODO, disable for the moment    TEST_THROW( needs.Get(aNeed,value), std::logic_error );
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
    needs.Set("foobar",trueValue);
    RCP<foobarClass> value;
    needs.Get("foobar",value);
    TEUCHOS_TEST_EQUALITY(trueValue,value, out, success);
  }

}//namespace <anonymous>
