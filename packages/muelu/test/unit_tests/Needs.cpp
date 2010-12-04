#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

namespace {

using Teuchos::RCP;
using Teuchos::rcp;
using MueLu::Needs;

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(NeedsObject, Constructor)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success
  out << "version: " << MueLu::Version() << std::endl;
  RCP<Needs> needs = rcp(new Needs() );
  TEUCHOS_TEST_EQUALITY(needs != Teuchos::null, true, out, success);
}

TEUCHOS_UNIT_TEST(NeedsObject, NeedRegistered)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  std::string aNeed = "knockNeed";
  TEUCHOS_TEST_EQUALITY(needs.IsRegistered(aNeed), false, out, success);
  needs.Want(aNeed);
  TEUCHOS_TEST_EQUALITY(needs.IsRegistered(aNeed), true, out, success);
}

TEUCHOS_UNIT_TEST(NeedsObject, GetCount_Exception)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  TEST_THROW( needs.GetCount("nonExistentNeed"), std::logic_error );
}

TEUCHOS_UNIT_TEST(NeedsObject, GetCount)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  std::string aNeed = "knockNeed";
  needs.Want(aNeed);
  needs.Want(aNeed);
  TEUCHOS_TEST_EQUALITY(needs.GetCount(aNeed), 2, out, success);
}

TEUCHOS_UNIT_TEST(NeedsObject, Get_Exception)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  double value=0;
  TEST_THROW( needs.Get("nonExistentNeed",value), std::logic_error );
}

TEUCHOS_UNIT_TEST(NeedsObject, SaveAndGet)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  std::string aNeed = "knockNeed";
  double trueValue = 42;
  needs.Save(aNeed,trueValue);
  double expectedValue = 0;
  needs.Get(aNeed,expectedValue);
  TEUCHOS_TEST_EQUALITY(trueValue,expectedValue, out, success);
}

TEUCHOS_UNIT_TEST(NeedsObject, Release_Exception)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  double value;
  TEST_THROW( needs.Release("nonExistentNeed",value), std::logic_error );
}

TEUCHOS_UNIT_TEST(NeedsObject, Release)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  std::string aNeed = "knockNeed";
  double trueValue = 42;
  needs.Save(aNeed,trueValue);
  needs.Want(aNeed);
  needs.Want(aNeed);
  double value = 0;
  needs.Release(aNeed,value);
  TEUCHOS_TEST_EQUALITY(trueValue,value, out, success);
  TEUCHOS_TEST_EQUALITY(needs.GetCount(aNeed),1, out, success);
  value = 0;
  needs.Release(aNeed,value);
  //try to get the need one too many times
  TEST_THROW( needs.Get(aNeed,value), std::logic_error );
}

class foobarClass {
  public:
  foobarClass() {}
  virtual ~foobarClass() {}
};

TEUCHOS_UNIT_TEST(NeedsObject, nonPOD)
{
  out << "version: " << MueLu::Version() << std::endl;
  Needs needs = Needs();
  RCP<foobarClass> trueValue = rcp(new foobarClass);
  needs.Save("foobar",trueValue);
  RCP<foobarClass> value;
  needs.Get("foobar",value);
  TEUCHOS_TEST_EQUALITY(trueValue,value, out, success);
}

}//namespace <anonymous>
