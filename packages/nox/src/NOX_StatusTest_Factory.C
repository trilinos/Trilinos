// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
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
buildStatusTests(const std::string& file_name , const NOX::Utils& u,
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
	       std::map<std::string, Teuchos::RCP<NOX::StatusTest::Generic> >* 
	       tagged_tests) const
{ 

  int number_of_tests = get<int>(p, "Number of Tests");
  
  std::string combo_type_string = get<std::string>(p, "Combo Type");
  NOX::StatusTest::Combo::ComboType combo_type;
  if (combo_type_string == "AND")
    combo_type = NOX::StatusTest::Combo::AND;
  else if (combo_type_string == "OR")
    combo_type = NOX::StatusTest::Combo::OR;
  else{
    std::string msg = 
      "Error - The \"Combo Type\" must be \"AND\" or \"OR\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
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
  double tolerance = p.get("Tolerance", 1.0e-8);
  
  // Norm Type
  std::string norm_type_string = p.get("Norm Type", "Two Norm");
  NOX::Abstract::Vector::NormType norm_type = NOX::Abstract::Vector::TwoNorm;
  if (norm_type_string == "Two Norm")
    norm_type = NOX::Abstract::Vector::TwoNorm;
  else if (norm_type_string == "One Norm")
    norm_type = NOX::Abstract::Vector::OneNorm;
  else if (norm_type_string == "Max Norm")
    norm_type = NOX::Abstract::Vector::MaxNorm;
  else {
    std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
  // Scale Type
  std::string scale_type_string = p.get("Scale Type", "Unscaled");
  NOX::StatusTest::NormF::ScaleType scale_type = 
    NOX::StatusTest::NormF::Unscaled;
  if (scale_type_string == "Unscaled")
    scale_type = NOX::StatusTest::NormF::Unscaled;
  else if (scale_type_string == "Scaled")
    scale_type = NOX::StatusTest::NormF::Scaled;
  else {
    std::string msg = "\"Scale Type\" must be either \"Unscaled\" or \"Scaled\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  // Relative or absoltue tolerance (relative requires f_0)
  bool use_relative_tolerance = false;
  Teuchos::RCP<NOX::Abstract::Group> group;
  if (isParameterType< RCP<NOX::Abstract::Group> >(p, "Initial Guess")) {
    group = get< RCP<NOX::Abstract::Group> >(p, "Initial Guess");
    use_relative_tolerance = true;
  }

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
buildNormUpdateTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  double tolerance = p.get("Tolerance", 1.0e-3);
  
  // Norm Type
  std::string norm_type_string = p.get("Norm Type", "Two Norm");
  NOX::Abstract::Vector::NormType norm_type = NOX::Abstract::Vector::TwoNorm;
  if (norm_type_string == "Two Norm")
    norm_type = NOX::Abstract::Vector::TwoNorm;
  else if (norm_type_string == "One Norm")
    norm_type = NOX::Abstract::Vector::OneNorm;
  else if (norm_type_string == "Max Norm")
    norm_type = NOX::Abstract::Vector::MaxNorm;
  else {
    std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
  // Scale Type
  std::string scale_type_string = p.get("Scale Type", "Unscaled");
  NOX::StatusTest::NormUpdate::ScaleType scale_type = 
    NOX::StatusTest::NormUpdate::Unscaled;
  if (scale_type_string == "Unscaled")
    scale_type = NOX::StatusTest::NormUpdate::Unscaled;
  else if (scale_type_string == "Scaled")
    scale_type = NOX::StatusTest::NormUpdate::Scaled;
  else {
    std::string msg = "\"Scale Type\" must be either \"Unscaled\" or \"Scaled\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  Teuchos::RCP<NOX::StatusTest::NormUpdate> status_test =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(tolerance, norm_type, 
						 scale_type));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildNormWRMSTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  double bdf_multiplier = p.get("BDF Multiplier", 1.0);
  double tolerance = p.get("Tolerance", 1.0);
  double alpha = p.get("Alpha", 1.0);
  double beta = p.get("Beta", 0.5);
  double rel_tol = p.get("Relative Tolerance", 1.0e-5);
  bool disable_implicit_weighting = p.get("Disable Implicit Weighting", true);

  bool abs_tol_is_vector = false;
  Teuchos::RCP<const NOX::Abstract::Vector> abs_tol_vector;
  double abs_tol = 1.0;
  if (isParameterType< RCP<const NOX::Abstract::Vector> >
      (p, "Absolute Tolerance")) {
    abs_tol_vector = get< Teuchos::RCP<const NOX::Abstract::Vector> >
      (p, "Absolute Tolerance");
    abs_tol_is_vector = true;
  }
  else {
    abs_tol = p.get("Absolute Tolerance", 1.0e-8);
  }
    
  RCP<NOX::StatusTest::NormWRMS> status_test;

  if (abs_tol_is_vector)
    status_test = rcp(new NOX::StatusTest::NormWRMS(rel_tol,
						    abs_tol_vector,
						    bdf_multiplier,
						    tolerance,
						    alpha,
						    beta,
						    disable_implicit_weighting)
		      );
  else
    status_test = rcp(new NOX::StatusTest::NormWRMS(rel_tol,
						    abs_tol,
						    bdf_multiplier,
						    tolerance,
						    alpha,
						    beta,
						    disable_implicit_weighting)
		      );

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildFiniteValueTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  std::string vector_type_string = p.get("Vector Type","F Vector");
  std::string norm_type_string = p.get("Norm Type", "Two Norm");

  NOX::StatusTest::FiniteValue::VectorType vector_type = 
    NOX::StatusTest::FiniteValue::FVector;
  NOX::Abstract::Vector::NormType norm_type = NOX::Abstract::Vector::TwoNorm;

  if (vector_type_string == "F Vector")
    vector_type = NOX::StatusTest::FiniteValue::FVector;
  else if (vector_type_string == "Solution Vector")
    vector_type = NOX::StatusTest::FiniteValue::SolutionVector;
  else {
    std::string msg = "\"Vector Type\" must be either \"F Vector\" or \"Solution Vector\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
  if (norm_type_string == "Two Norm")
    norm_type = NOX::Abstract::Vector::TwoNorm;
  else if (vector_type_string == "One Norm")
    norm_type = NOX::Abstract::Vector::OneNorm;
  else if (vector_type_string == "Max Norm")
    norm_type = NOX::Abstract::Vector::MaxNorm;
  else {
    std::string msg = "\"Norm Type\" must be either \"Two Norm\", \"One Norm\", or \"Max Norm\"!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
  RCP<NOX::StatusTest::FiniteValue> status_test = 
    rcp(new NOX::StatusTest::FiniteValue(vector_type, norm_type));

  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildDivergenceTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  double tolerance = p.get("Tolerance", 1.0e+12);
  int iterations = p.get("Consecutive Iterations", 1);
  
  RCP<NOX::StatusTest::Divergence> status_test = 
    rcp(new NOX::StatusTest::Divergence(tolerance, iterations));
  
  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildStagnationTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  double tolerance = p.get("Tolerance", 1.0e+12);
  int iterations = p.get("Consecutive Iterations", 1);
  
  RCP<NOX::StatusTest::Stagnation> status_test = 
    rcp(new NOX::StatusTest::Stagnation(iterations, tolerance));
  
  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildMaxItersTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  int max_iters = get<int>(p, "Maximum Iterations");
  
  RCP<NOX::StatusTest::MaxIters> status_test = 
    rcp(new NOX::StatusTest::MaxIters(max_iters, &u));
  
  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildRelativeNormFTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  double tolerance = p.get("Tolerance", 1.0e-8);
  bool scale_by_length = p.get("Scale Norms by Length", false);

  RCP<NOX::StatusTest::RelativeNormF> status_test;

  status_test = rcp(new NOX::StatusTest::RelativeNormF(tolerance, 
						       scale_by_length,
						       &u));
  
  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildNStepTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  int num_iters = p.get<int>("Number of Nonlinear Iterations", 1);
  int num_ramping_steps = p.get<int>("Number of Initial Ramping Steps", 0);
  int num_ramping_iters = p.get<int>("Number of Nonlinear Iterations in Ramping Phase", 10);

  RCP<NOX::StatusTest::NStep> status_test;

  status_test = rcp(new NOX::StatusTest::NStep(num_iters, num_ramping_steps,
					       num_ramping_iters));
  
  return status_test;
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::StatusTest::Generic> NOX::StatusTest::Factory::
buildUserDefinedTest(Teuchos::ParameterList& p, const NOX::Utils& u) const
{
  RCP<NOX::StatusTest::Generic> status_test;
  
  if (isParameterType< RCP<NOX::StatusTest::Generic> >(p, "User Status Test"))
    status_test = get< RCP<NOX::StatusTest::Generic> >(p, "User Status Test");
  else {
    std::string msg = "Error - NOX::StatusTest::Factory::buildUserDefinedTest() - a user defined status test has been selected, but the test has not been supplied as an RCP<NOX::StatusTest::Generic> in the parameter list.  please make sure it is set as a \"Generic\" object in the parameter list.";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
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
  if ( (isParameterType<std::string>(p, "Tag")) && (tagged_tests != NULL) ) {
    (*tagged_tests)[getParameter<std::string>(p, "Tag")] = test;
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
