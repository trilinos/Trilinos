#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

namespace MueLuTests {

  using MueLu::Needs;

  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST(Needs, Constructor)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Needs> needs = rcp(new Needs() );
    TEST_EQUALITY(needs != Teuchos::null, true);
  }

  TEUCHOS_UNIT_TEST(Needs, NeedRequested)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    TEST_EQUALITY(needs.IsRequested(aNeed,NULL), false);
    needs.Request(aNeed,NULL);
    TEST_EQUALITY(needs.IsRequested(aNeed,NULL), true);
  }

  TEUCHOS_UNIT_TEST(Needs, ValueIsAvailable)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    TEST_EQUALITY(needs.IsAvailable(aNeed,NULL), false);
    needs.SetData(aNeed,42);
    TEST_EQUALITY(needs.IsAvailable(aNeed,NULL), true);
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
    needs.Request(aNeed,NULL);
    needs.Request(aNeed,NULL);
    TEST_EQUALITY(needs.NumRequests(aNeed), 2);
  }

  TEUCHOS_UNIT_TEST(Needs, Get_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    double value=0;
    TEST_THROW( needs.GetData("nonExistentNeed",value), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, SetAndGet)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.SetData(aNeed,trueValue);
    double expectedValue = 0;
    needs.GetData(aNeed,expectedValue);
    TEST_EQUALITY(trueValue,expectedValue);
  }

  TEUCHOS_UNIT_TEST(Needs, Release_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    TEST_THROW( needs.Release("nonExistentNeed",NULL), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, Release_Without_Request)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.SetData(aNeed,trueValue);
    //    double expectedValue = 0;
    //JG TODO
//     TEST_THROW( needs.Get(aNeed,expectedValue), MueLu::Exceptions::RuntimeError );
//     TEST_THROW( needs.Release(aNeed), MueLu::Exceptions::RuntimeError );
  }

  TEUCHOS_UNIT_TEST(Needs, Release)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.SetData(aNeed,trueValue);
    needs.Request(aNeed,NULL);         // TODO: write new test
    needs.Request(aNeed,NULL);
    double value = 0;
    needs.GetData(aNeed,value);
    needs.Release(aNeed,NULL);
    TEST_EQUALITY(trueValue,value);
    TEST_EQUALITY(needs.NumRequests(aNeed),1);
    value = 0;
    needs.GetData(aNeed,value);
    needs.Release(aNeed,NULL);
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
    needs.SetData("foobar",trueValue);
    RCP<foobarClass> value;
    needs.GetData("foobar",value);
    TEST_EQUALITY(trueValue,value);
  }

}//namespace MueLuTests
