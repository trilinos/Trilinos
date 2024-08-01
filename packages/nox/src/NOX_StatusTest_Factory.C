// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListHelpers.hpp"
#endif

#include "NOX_Common.H"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Factory.H"
#include "NOX_Utils.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"

// Concrete StatusTest Objects
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_StatusTest_Divergence.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_RelativeNormF.H"
#include "NOX_StatusTest_NStep.H"

#include <sstream>

using namespace Teuchos;

// ************************************************************************
// ************************************************************************
NOX::StatusTest::Factory::Factory()
{ }

// ************************************************************************
// ************************************************************************
NOX::StatusTest::Factory::~Factory()
{ }

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildStatusTests(const std::string& /* file_name */ , const NOX::Utils& u,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
                 tagged_tests) const
{
  Teuchos::RCP<NOX::StatusTest::Generic> status_tests;

#ifdef HAVE_TEUCHOS_EXTENDED
  Teuchos::ParameterList param_list;
  Teuchos::updateParametersFromXmlFile("input.xml", ptrFromRef(param_list));
  status_tests = this->buildStatusTests(param_list, u, tagged_tests);
#else
  std::string msg = "Error - Teuchos Extended Support must be enabled to use the xml reader for parameter lists.  Please rebuild the Trilinos Teuchos library with exteded support enabled.";
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
#endif

  return status_tests;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildStatusTests(Teuchos::ParameterList& p, const NOX::Utils& u,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
                 tagged_tests) const
{
  Teuchos::RCP<NOX::StatusTest::Generic> status_test;

  std::string test_type = "???";

  if (isParameterType<std::string>(p, "Test Type"))
    test_type = get<std::string>(p, "Test Type");
  else {
    std::string msg = "Error - The \"Test Type\" is a required parameter in the NOX::StatusTest::Factory!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  if (test_type == "Combo")
    status_test = this->buildComboTest(p, u, tagged_tests);
  else if (test_type == "NormF")
    status_test = this->buildNormFTest(p, u);
  else if (test_type == "NormUpdate")
    status_test = this->buildNormUpdateTest(p, u);
  else if (test_type == "NormWRMS")
    status_test = this->buildNormWRMSTest(p, u);
  else if (test_type == "FiniteValue")
    status_test = this->buildFiniteValueTest(p, u);
  else if (test_type == "MaxIters")
    status_test = this->buildMaxItersTest(p, u);
  else if (test_type == "Divergence")
    status_test = this->buildDivergenceTest(p, u);
  else if (test_type == "Stagnation")
    status_test = this->buildStagnationTest(p, u);
  else if (test_type == "RelativeNormF")
    status_test = this->buildRelativeNormFTest(p, u);
  else if (test_type == "NStep")
    status_test = this->buildNStepTest(p, u);
  else if (test_type == "User Defined")
    status_test = this->buildUserDefinedTest(p, u);
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
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildComboTest(Teuchos::ParameterList& p, const NOX::Utils& u,
               std::map<std::string,Teuchos::RCP<NOX::StatusTest::Generic>>* tagged_tests) const
{
  // Validation is tricky due to number of tests being runtime
  // parameter.
  int number_of_tests = p.get<int>("Number of Tests"); // Must be set by user

  Teuchos::ParameterList validParams;
  validParams.set("Test Type","Combo");
  validParams.set<int>("Number of Tests",number_of_tests);
  for (int t=0; t < number_of_tests; ++t) {
    std::ostringstream subtest_name;
    subtest_name << "Test " << t;
    validParams.sublist(subtest_name.str());
  }
  // ROGER disable recursive
  Teuchos::setStringToIntegralParameter<NOX::StatusTest::Combo::ComboType>
    ("Combo Type", "AND", "Type of combination test to use, \"AND\" or \"OR\".",
     Teuchos::tuple<std::string> ("AND","OR"),
     Teuchos::tuple<NOX::StatusTest::Combo::ComboType>(NOX::StatusTest::Combo::AND,
                                                       NOX::StatusTest::Combo::OR),
     &validParams);
  validParams.set("Tag","");
  int validation_depth = 0; // Do not validate the sublists, they are
                            // validated separately for each test.
  p.validateParametersAndSetDefaults(validParams,validation_depth);

  auto combo_type = Teuchos::getIntegralValue<NOX::StatusTest::Combo::ComboType>(p,"Combo Type");

  RCP<NOX::StatusTest::Combo> combo_test =
    rcp(new NOX::StatusTest::Combo(combo_type, &u));

  for (int i=0; i < number_of_tests; ++i) {
    std::ostringstream subtest_name;
    subtest_name << "Test " << i;
    ParameterList& subtest_list = p.sublist(subtest_name.str(), true);

    RCP<NOX::StatusTest::Generic> subtest =
      this->buildStatusTests(subtest_list, u, tagged_tests);

    combo_test->addStatusTest(subtest);
  }

  return combo_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildNormFTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","NormF");
  validParams.set<double>("Tolerance",1.0e-8);
  Teuchos::setStringToIntegralParameter<NOX::Abstract::Vector::NormType>
    ("Norm Type", "Two Norm", "Type of norm to use.",
     Teuchos::tuple<std::string>("Two Norm","One Norm","Max Norm"),
     Teuchos::tuple<NOX::Abstract::Vector::NormType>(NOX::Abstract::Vector::TwoNorm,
                                                     NOX::Abstract::Vector::OneNorm,
                                                     NOX::Abstract::Vector::MaxNorm),
     &validParams);
  Teuchos::setStringToIntegralParameter<NOX::StatusTest::NormF::ScaleType>
    ("Scale Type", "Unscaled", "Whether to use scaled or unscaled norms.",
     Teuchos::tuple<std::string> ("Unscaled","Scaled"),
     Teuchos::tuple<NOX::StatusTest::NormF::ScaleType>(NOX::StatusTest::NormF::Unscaled,
                                                       NOX::StatusTest::NormF::Scaled),
     &validParams);
  validParams.set<Teuchos::RCP<NOX::Abstract::Group>>("Initial Guess",Teuchos::null);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  double tolerance = p.get<double>("Tolerance");
  auto norm_type = Teuchos::getIntegralValue<NOX::Abstract::Vector::NormType>(p,"Norm Type");
  auto scale_type = Teuchos::getIntegralValue<NOX::StatusTest::NormF::ScaleType>(p,"Scale Type");

  // Relative or absoltue tolerance (relative requires f_0)
  bool use_relative_tolerance = false;
  auto group = p.get<Teuchos::RCP<NOX::Abstract::Group>>("Initial Guess");
  if (nonnull(group))
    use_relative_tolerance = true;

  RCP<NOX::StatusTest::NormF> status_test;

  if (use_relative_tolerance)
    status_test = rcp(new NOX::StatusTest::NormF(*group,
                                                 tolerance,
                                                 norm_type,
                                                 scale_type,
                                                 &u));
  else
    status_test = rcp(new NOX::StatusTest::NormF(tolerance,
                                                 norm_type,
                                                 scale_type,
                                                 &u));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildNormUpdateTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","NormUpdate");
  validParams.set<double>("Tolerance",1.0e-3);
  Teuchos::setStringToIntegralParameter<NOX::Abstract::Vector::NormType>
    ("Norm Type", "Two Norm", "Type of norm to use.",
     Teuchos::tuple<std::string>("Two Norm","One Norm","Max Norm"),
     Teuchos::tuple<NOX::Abstract::Vector::NormType>(NOX::Abstract::Vector::TwoNorm,
                                                     NOX::Abstract::Vector::OneNorm,
                                                     NOX::Abstract::Vector::MaxNorm),
     &validParams);
  Teuchos::setStringToIntegralParameter<NOX::StatusTest::NormUpdate::ScaleType>
    ("Scale Type", "Unscaled", "Whether to scaled or unscaled norms.",
     Teuchos::tuple<std::string> ("Unscaled","Scaled"),
     Teuchos::tuple<NOX::StatusTest::NormUpdate::ScaleType>(NOX::StatusTest::NormUpdate::Unscaled,
                                                            NOX::StatusTest::NormUpdate::Scaled),
     &validParams);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  double tolerance = p.get<double>("Tolerance");
  auto norm_type = Teuchos::getIntegralValue<NOX::Abstract::Vector::NormType>(p,"Norm Type");
  auto scale_type = Teuchos::getIntegralValue<NOX::StatusTest::NormUpdate::ScaleType>(p,"Scale Type");

  Teuchos::RCP<NOX::StatusTest::NormUpdate> status_test =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(tolerance,norm_type,scale_type));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildNormWRMSTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  // PList validation is tricky because the "Absolute Tolerance"
  // parameter can be a scalar value or a vector value. Have to check
  // the type up front to build the correct validator list.
  bool abs_tol_is_vector = false;
  if (Teuchos::isParameterType<Teuchos::RCP<const NOX::Abstract::Vector>>(p, "Absolute Tolerance"))
    abs_tol_is_vector = true;

  Teuchos::ParameterList validParams;
  validParams.set("Test Type","NormWRMS");
  validParams.set<double>("BDF Multiplier",1.0);
  validParams.set<double>("Tolerance",1.0);
  validParams.set<double>("Alpha",1.0);
  validParams.set<double>("Beta",0.5);
  validParams.set<double>("Relative Tolerance",1.0e-5);
  validParams.set<bool>("Disable Implicit Weighting",true);
  if (abs_tol_is_vector)
    validParams.set<Teuchos::RCP<const NOX::Abstract::Vector>>("Absolute Tolerance",Teuchos::null);
  else
    validParams.set<double>("Absolute Tolerance",1.0e-8);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  double bdf_multiplier = p.get<double>("BDF Multiplier");
  double tolerance = p.get<double>("Tolerance");
  double alpha = p.get<double>("Alpha");
  double beta = p.get<double>("Beta");
  double rel_tol = p.get<double>("Relative Tolerance");
  bool disable_implicit_weighting = p.get<bool>("Disable Implicit Weighting");

  Teuchos::RCP<const NOX::Abstract::Vector> abs_tol_vector;
  double abs_tol = 1.0;
  if (abs_tol_is_vector)
    abs_tol_vector = get< Teuchos::RCP<const NOX::Abstract::Vector> >(p, "Absolute Tolerance");
  else
    abs_tol = p.get<double>("Absolute Tolerance");

  RCP<NOX::StatusTest::NormWRMS> status_test;

  if (abs_tol_is_vector)
    status_test = rcp(new NOX::StatusTest::NormWRMS(rel_tol,
                                                    abs_tol_vector,
                                                    bdf_multiplier,
                                                    tolerance,
                                                    alpha,
                                                    beta,
                                                    disable_implicit_weighting));
  else
    status_test = rcp(new NOX::StatusTest::NormWRMS(rel_tol,
                                                    abs_tol,
                                                    bdf_multiplier,
                                                    tolerance,
                                                    alpha,
                                                    beta,
                                                    disable_implicit_weighting));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildFiniteValueTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","FiniteValue");
  Teuchos::setStringToIntegralParameter<NOX::StatusTest::FiniteValue::VectorType>
    ("Vector Type", "F Vector", "Which vector to check, F or X.",
     Teuchos::tuple<std::string>("F Vector","Solution Vector"),
     Teuchos::tuple<NOX::StatusTest::FiniteValue::VectorType>(NOX::StatusTest::FiniteValue::FVector,
                                                              NOX::StatusTest::FiniteValue::SolutionVector),
     &validParams);
  Teuchos::setStringToIntegralParameter<NOX::Abstract::Vector::NormType>
    ("Norm Type", "Two Norm", "Type of norm to use.",
     Teuchos::tuple<std::string>("Two Norm","One Norm","Max Norm"),
     Teuchos::tuple<NOX::Abstract::Vector::NormType>(NOX::Abstract::Vector::TwoNorm,
                                                     NOX::Abstract::Vector::OneNorm,
                                                     NOX::Abstract::Vector::MaxNorm),
     &validParams);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  auto vector_type = Teuchos::getIntegralValue<NOX::StatusTest::FiniteValue::VectorType>(p,"Vector Type");
  auto norm_type = Teuchos::getIntegralValue<NOX::Abstract::Vector::NormType>(p,"Norm Type");

  RCP<NOX::StatusTest::FiniteValue> status_test =
    rcp(new NOX::StatusTest::FiniteValue(vector_type, norm_type));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildDivergenceTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","Divergence");
  validParams.set<double>("Tolerance",1.0e+12);
  validParams.set<int>("Consecutive Iterations",1);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  double tolerance = p.get<double>("Tolerance");
  int iterations = p.get<int>("Consecutive Iterations");

  RCP<NOX::StatusTest::Divergence> status_test =
    rcp(new NOX::StatusTest::Divergence(tolerance, iterations));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildStagnationTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","Stagnation");
  validParams.set<double>("Tolerance",1.0e+12);
  validParams.set<int>("Consecutive Iterations",1);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  double tolerance = p.get<double>("Tolerance");
  int iterations = p.get<int>("Consecutive Iterations");

  RCP<NOX::StatusTest::Stagnation> status_test =
    rcp(new NOX::StatusTest::Stagnation(iterations, tolerance));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildMaxItersTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","MaxIters");
  validParams.set("Maximum Iterations",20);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  int max_iters = p.get<int>("Maximum Iterations");

  RCP<NOX::StatusTest::MaxIters> status_test =
    rcp(new NOX::StatusTest::MaxIters(max_iters, &u));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildRelativeNormFTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","RelativeNormF");
  validParams.set("Tolerance",1.0e-8);
  validParams.set("Scale Norms by Length",false);
  Teuchos::setStringToIntegralParameter<NOX::Abstract::Vector::NormType>
    ("Norm Type", "Two Norm", "Type of norm to use",
     Teuchos::tuple<std::string>("Two Norm","One Norm","Max Norm"),
     Teuchos::tuple<NOX::Abstract::Vector::NormType>(NOX::Abstract::Vector::TwoNorm,
                                                     NOX::Abstract::Vector::OneNorm,
                                                     NOX::Abstract::Vector::MaxNorm),
     &validParams);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  double tolerance = p.get<double>("Tolerance");
  bool scale_by_length = p.get<bool>("Scale Norms by Length");
  auto norm_type = Teuchos::getIntegralValue<NOX::Abstract::Vector::NormType>(p,"Norm Type");

  RCP<NOX::StatusTest::RelativeNormF> status_test;

  status_test = rcp(new NOX::StatusTest::RelativeNormF(tolerance,
                                                       scale_by_length,
                                                       &u,
                                                       norm_type));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildNStepTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","NStep");
  validParams.set("Number of Nonlinear Iterations",1);
  validParams.set("Number of Initial Ramping Steps",0);
  validParams.set("Number of Nonlinear Iterations in Ramping Phase",10);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  int num_iters = p.get<int>("Number of Nonlinear Iterations");
  int num_ramping_steps = p.get<int>("Number of Initial Ramping Steps");
  int num_ramping_iters = p.get<int>("Number of Nonlinear Iterations in Ramping Phase");

  RCP<NOX::StatusTest::NStep> status_test;

  status_test = rcp(new NOX::StatusTest::NStep(num_iters, num_ramping_steps,
                           num_ramping_iters));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildUserDefinedTest(Teuchos::ParameterList& p, const NOX::Utils& /* u */) const
{
  Teuchos::ParameterList validParams;
  validParams.set("Test Type","UserDefined");
  validParams.set<Teuchos::RCP<NOX::StatusTest::Generic>>("User Status Test",Teuchos::null);
  validParams.set("Tag","");
  p.validateParametersAndSetDefaults(validParams);

  auto status_test = p.get<Teuchos::RCP<NOX::StatusTest::Generic>>("User Status Test");

  return status_test;
}

// ************************************************************************
// ************************************************************************
bool NOX::StatusTest::Factory::
checkAndTagTest(const Teuchos::ParameterList& p,
                const Teuchos::RCP<NOX::StatusTest::Generic>& test,
                std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
                tagged_tests) const
{
  auto tag_name = p.get<std::string>("Tag");
  if ( (tag_name != "") && (tagged_tests != nullptr) ) {
    (*tagged_tests)[tag_name] = test;
    return true;
  }

  return false;
}

// ************************************************************************
// ************************************************************************
// Nonmember function
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::
buildStatusTests(const std::string& file_name, const NOX::Utils& utils,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
                 tagged_tests)
{
  NOX::StatusTest::Factory factory;
  return factory.buildStatusTests(file_name, utils, tagged_tests);
}

// ************************************************************************
// ************************************************************************
// Nonmember function
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::
buildStatusTests(Teuchos::ParameterList& p, const NOX::Utils& utils,
                 std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >*
                 tagged_tests)
{
  NOX::StatusTest::Factory factory;
  return factory.buildStatusTests(p, utils, tagged_tests);
}

// ************************************************************************
// ************************************************************************
