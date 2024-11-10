// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

#include "LOCA_StatusTest_Factory.H"

// Concrete StatusTest Objects
#include "LOCA_StatusTest_Combo.H"
#include "LOCA_StatusTest_MaxIters.H"

using namespace Teuchos;

// ************************************************************************
// ************************************************************************
LOCA::StatusTest::Factory::Factory()
{ }

// ************************************************************************
// ************************************************************************
LOCA::StatusTest::Factory::~Factory()
{ }

// ************************************************************************
// ************************************************************************
Teuchos::RCP<LOCA::StatusTest::Abstract> LOCA::StatusTest::Factory::
buildStatusTests(const std::string& file_name ,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
             std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests) const
{
  Teuchos::RCP<LOCA::StatusTest::Abstract> status_tests;

#ifdef HAVE_TEUCHOS_EXTENDED
  Teuchos::ParameterList param_list;
  const Teuchos::Ptr<Teuchos::ParameterList> param_list_ptr(&param_list);
  Teuchos::updateParametersFromXmlFile("input.xml", param_list_ptr);
  status_tests = this->buildStatusTests(param_list, globalData, tagged_tests);
#else
  std::string msg = "Error - Teuchos Extended Support must be enabled to use the xml reader for parameter lists.  Please rebuild the Trilinos Teuchos library with --enable-teuchos-extended in teh configure script.";
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
#endif

  return status_tests;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<LOCA::StatusTest::Abstract> LOCA::StatusTest::Factory::
buildStatusTests(Teuchos::ParameterList& p,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
          std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests) const
{
  Teuchos::RCP<LOCA::StatusTest::Abstract> status_test;

  std::string test_type = "???";

  if (isParameterType<std::string>(p, "Test Type"))
    test_type = get<std::string>(p, "Test Type");
  else {
    std::string msg = "Error - The \"Test Type\" is a required parameter in the LOCA::StatusTest::Factory!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  if (test_type == "Combo")
    status_test = this->buildComboTest(p, globalData, tagged_tests);
  else if (test_type == "MaxIters")
    status_test = this->buildMaxItersTest(p, globalData);
  else {
    std::ostringstream msg;
    msg << "Error - the test type \"" << test_type << "\" is invalid!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  this->checkAndTagTest(p, status_test, tagged_tests);

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<LOCA::StatusTest::Abstract> LOCA::StatusTest::Factory::
buildComboTest(Teuchos::ParameterList& p,
               const Teuchos::RCP<const LOCA::GlobalData> & globalData,
           std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
           tagged_tests) const
{

  int number_of_tests = get<int>(p, "Number of Tests");

  std::string combo_type_string = get<std::string>(p, "Combo Type");
  LOCA::StatusTest::Combo::ComboType combo_type;
  if (combo_type_string == "AND")
    combo_type = LOCA::StatusTest::Combo::AND;
  else if (combo_type_string == "OR")
    combo_type = LOCA::StatusTest::Combo::OR;
  else{
    std::string msg =
      "Error - The \"Combo Type\" must be \"AND\" or \"OR\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  RCP<LOCA::StatusTest::Combo> combo_test =
    rcp(new LOCA::StatusTest::Combo(combo_type, globalData));

  for (int i=0; i < number_of_tests; ++i) {
    ostringstream subtest_name;
    subtest_name << "Test " << i;
    ParameterList& subtest_list = p.sublist(subtest_name.str(), true);

    RCP<LOCA::StatusTest::Abstract> subtest =
      this->buildStatusTests(subtest_list, globalData, tagged_tests);

    combo_test->addStatusTest(subtest);
  }

  return combo_test;
}
// ************************************************************************
// ************************************************************************
Teuchos::RCP<LOCA::StatusTest::Abstract> LOCA::StatusTest::Factory::
buildMaxItersTest(Teuchos::ParameterList& p,
                  const Teuchos::RCP<const LOCA::GlobalData> & globalData) const
{
  int max_iters = get<int>(p, "Maximum Iterations");
  bool return_failed_on_max_steps = p.get<bool>("Return failed on max steps", true);

  RCP<LOCA::StatusTest::MaxIters> status_test =
    rcp(new LOCA::StatusTest::MaxIters(max_iters, return_failed_on_max_steps, globalData));

  return status_test;
}

// ************************************************************************
// ************************************************************************
bool LOCA::StatusTest::Factory::
checkAndTagTest(const Teuchos::ParameterList& p,
        const Teuchos::RCP<LOCA::StatusTest::Abstract>& test,
        std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
        tagged_tests) const
{
  if ( (isParameterType<std::string>(p, "Tag")) && (tagged_tests != NULL) ) {
    (*tagged_tests)[getParameter<std::string>(p, "Tag")] = test;
    return true;
  }

  return false;
}

// ************************************************************************
// ************************************************************************
// Nonmember function
Teuchos::RCP<LOCA::StatusTest::Abstract> LOCA::StatusTest::
buildStatusTests(const std::string& file_name,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
             std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests)
{
  LOCA::StatusTest::Factory factory;
  return factory.buildStatusTests(file_name, globalData, tagged_tests);
}

// ************************************************************************
// ************************************************************************
// Nonmember function
Teuchos::RCP<LOCA::StatusTest::Abstract> LOCA::StatusTest::
buildStatusTests(Teuchos::ParameterList& p,
                 const Teuchos::RCP<const LOCA::GlobalData> & globalData,
         std::map<std::string, Teuchos::RCP<LOCA::StatusTest::Abstract> >*
         tagged_tests)
{
  LOCA::StatusTest::Factory factory;
  return factory.buildStatusTests(p, globalData, tagged_tests);
}

// ************************************************************************
// ************************************************************************
