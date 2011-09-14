#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

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
    needs.Request(aNeed,NULL);
    needs.SetData(aNeed,42,NULL);
    TEST_EQUALITY(needs.IsAvailable(aNeed,NULL), true);
  }

  TEUCHOS_UNIT_TEST(Needs, NumRequests_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    TEST_THROW( needs.NumRequests("nonExistentNeed",NULL), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, NumRequests)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    needs.Request(aNeed,NULL);
    needs.Request(aNeed,NULL);
    TEST_EQUALITY(needs.NumRequests(aNeed,NULL), 2);
  }

  TEUCHOS_UNIT_TEST(Needs, Get_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    double value=0;
    TEST_THROW( needs.GetData("nonExistentNeed",value,NULL), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, SetAndGet)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Request(aNeed,NULL);
    needs.SetData(aNeed,trueValue,NULL);
    double expectedValue = 0;
    needs.GetData(aNeed,expectedValue,NULL);
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
    needs.SetData(aNeed,trueValue,NULL);
    double expectedValue = 0;
    TEST_THROW( expectedValue = needs.GetData<double>(aNeed,NULL), MueLu::Exceptions::RuntimeError );
    TEST_THROW( needs.Release(aNeed,NULL), MueLu::Exceptions::RuntimeError );

	//    RCP<MueLu::TentativePFactory> fac = rcp(new MueLu::TentativePFactory() );
	//    needs.SetData<double>("test", trueValue, &fac);
	//    TEST_THROW( expectedValue = needs.GetData<double>("test",&fac), MueLu::Exceptions::RuntimeError );
	//    TEST_THROW( needs.Release("test",&fac), MueLu::Exceptions::RuntimeError );

  }

  TEUCHOS_UNIT_TEST(Needs, Release)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs = Needs();
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Request(aNeed,NULL);         // TODO: write new test
    needs.Request(aNeed,NULL);
    needs.SetData(aNeed,trueValue,NULL);
    double value = 0;
    needs.GetData(aNeed,value,NULL);
    needs.Release(aNeed,NULL);
    TEST_EQUALITY(trueValue,value);
    TEST_EQUALITY(needs.NumRequests(aNeed,NULL),1);
    value = 0;
    needs.GetData(aNeed,value,NULL);
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
    needs.Request("foobar",NULL);
    needs.SetData("foobar",trueValue,NULL);
    RCP<foobarClass> value;
    needs.GetData("foobar",value,NULL);
    TEST_EQUALITY(trueValue,value);
  }

}//namespace MueLuTests
