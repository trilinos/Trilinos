// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_util/stk_config.h>

#ifdef STK_HAVE_BOOST

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <sys/types.h>                  // for int64_t
#include <complex>                      // for complex, operator<<, etc
#include <exception>                    // for exception
#include <iostream>                     // for basic_ostream::operator<<, etc
#include <map>                          // for operator==, etc
#include <stk_util/util/ParameterList.hpp>  // for ParameterList, Type, etc
#include <string>                       // for string, basic_string
#include <vector>                       // for vector
#include "boost/any.hpp"                // for any_cast

namespace
{
  TEST(StkUtilTestForDocumentation, Parameters)
  {
    //BEGINParametersInit
    //+ INITIALIZATION
    std::vector<std::string> expected_name;
    std::vector<stk::util::ParameterType::Type> expected_type;
    
    //+ Scalar values of type double, float, int, int64_t, and string
    double pi = 3.14159;
    float e = 2.71828;
    int answer = 42;
    int64_t big_answer = 42000000000001;
    std::string team_name = "STK Transition Team";

    expected_name.push_back("PI");
    expected_type.push_back(stk::util::ParameterType::DOUBLE);
    expected_name.push_back("E");
    expected_type.push_back(stk::util::ParameterType::FLOAT);
    expected_name.push_back("Answer");
    expected_type.push_back(stk::util::ParameterType::INTEGER);
    expected_name.push_back("Answer_64");
    expected_type.push_back(stk::util::ParameterType::INT64);
    expected_name.push_back("TeamName");
    expected_type.push_back(stk::util::ParameterType::STRING);

    //+ vector of doubles
    std::vector<double> my_double_vector;
    my_double_vector.push_back(2.78); my_double_vector.push_back(5.30);
    my_double_vector.push_back(6.21);
    expected_name.push_back("some_doubles");
    expected_type.push_back(stk::util::ParameterType::DOUBLEVECTOR);
    
    //+ vector of floats
    std::vector<float> my_float_vector;
    my_float_vector.push_back(194.0); my_float_vector.push_back(-194.0);
    my_float_vector.push_back(47.0);  my_float_vector.push_back(92.0);
    expected_name.push_back("some_floats");
    expected_type.push_back(stk::util::ParameterType::FLOATVECTOR);
    
    //+ vector of ints
    std::vector<int> ages;
    ages.push_back(55); ages.push_back(49); ages.push_back(21); ages.push_back(19);
    expected_name.push_back("Ages");
    expected_type.push_back(stk::util::ParameterType::INTEGERVECTOR);
    
    //+ vector of int64_ts
    std::vector<int64_t> ages_64;
    ages_64.push_back(55); ages_64.push_back(49); ages_64.push_back(21); ages_64.push_back(19);
    expected_name.push_back("Ages_64");
    expected_type.push_back(stk::util::ParameterType::INT64VECTOR);
    
    //+ vector of strings
    std::vector<std::string> names;
    names.push_back("greg"); names.push_back("chloe"); names.push_back("tuffy");
    names.push_back("liberty"); names.push_back("I have spaces");
    expected_name.push_back("Names");
    expected_type.push_back(stk::util::ParameterType::STRINGVECTOR);
    //ENDParametersInit
    
    //BEGINParametersDefine
    //+ Define parameters...
    stk::util::ParameterList params;
    params.set_param("PI",           pi);
    params.set_param("E",            e);
    params.set_param("Answer",       answer);
    params.set_param("Answer_64",    big_answer);
    params.set_param("TeamName",     team_name);
    params.set_param("some_doubles", my_double_vector);
    params.set_param("some_floats",  my_float_vector); 
    params.set_param("Ages",         ages); 
    params.set_param("Ages_64",      ages_64); 
    params.set_param("Names",        names); 
    //ENDParametersDefine

    //BEGINParametersAccessingValues
    //+ Write parameters to stdout...
    params.write_parameter_list(std::cout);
      
    //+ Access parameters by name...
    size_t num_param = expected_name.size();
    for (size_t i=0; i < num_param; i++) {
      stk::util::Parameter &param = params.get_param(expected_name[i]);
      EXPECT_EQ(param.type, expected_type[i]);
    }

    //+ Extract some parameter values if know type:
    std::vector<int> pages = params.get_value<std::vector<int> >("Ages");
    for (size_t i=0; i < pages.size(); i++) {
      EXPECT_EQ(pages[i], ages[i]);
    }

    double my_pi = params.get_value<double>("PI");
    EXPECT_EQ(my_pi, pi);

    //+ Change value of an existing parameter
    params.set_value("Answer", 21);

    int new_answer = params.get_value<int>("Answer");
    EXPECT_EQ(new_answer, 21);
    
    {
      //+ Access a variable of unknown type...
      //+ The parameter uses boost::any to store the actual value.
      stk::util::Parameter &param = params.get_param("Answer");
      double value_as_double = 0.0;
      switch (param.type) {
      case stk::util::ParameterType::DOUBLE:
	value_as_double = boost::any_cast<double>(param.value);
	break;
      case stk::util::ParameterType::FLOAT:
	value_as_double = static_cast<double>(boost::any_cast<float>(param.value));
	break;
      case stk::util::ParameterType::INTEGER:
	value_as_double = static_cast<double>(boost::any_cast<int>(param.value));
	break;
      case stk::util::ParameterType::INT64:
	value_as_double = static_cast<double>(boost::any_cast<int64_t>(param.value));
	break;
      default:
	std::cerr << "ERROR: I can not convert 'Answers' value to a double\n";
	break;
      }
      EXPECT_EQ(static_cast<double>(new_answer), value_as_double);
    }

    {
      //+ Access a variable of unknown type without using boost::any_cast
      stk::util::Parameter &param = params.get_param("Answer");
      double value_as_double = 0.0;
      switch (param.type) {
      case stk::util::ParameterType::DOUBLE:
	value_as_double = params.get_value<double>("Answer");
	break;
      case stk::util::ParameterType::FLOAT:
	value_as_double = static_cast<double>(params.get_value<float>("Answer"));
	break;
      case stk::util::ParameterType::INTEGER:
	value_as_double = static_cast<double>(params.get_value<int>("Answer"));
	break;
      case stk::util::ParameterType::INT64:
	value_as_double = static_cast<double>(params.get_value<int64_t>("Answer"));
	break;
      default:
	std::cerr << "ERROR: I can not convert 'Answers' value to a double\n";
	break;
      }
      EXPECT_EQ(static_cast<double>(new_answer), value_as_double);
    }
    
    //ENDParametersAccessingValues

    //BEGINParametersErrors
    //+ If the requested parameter does not exist, 
    //+ an error message is printed to stderr and an invalid
    //+ parameter object is returned
    stk::util::Parameter no_exist = params.get_param("DoesNotExist");
    EXPECT_EQ(stk::util::ParameterType::INVALID, no_exist.type);
    
    //+ In this method of requesting a parameter, no error
    //+ message is printed if the parameter doesn't exist and
    //+ instead the returned iterator is equal to the end of the
    //+ parameter list.
    stk::util::ParameterMapType::iterator it = params.find("DoesNotExist");
    EXPECT_TRUE(it == params.end());
    
    //+ If the value of a non-existant parameter is requested,
    //+ an error message is printed and the value 0 is returned.
    double invalid_value = params.get_value<double>("DoesNotExist");
    EXPECT_EQ(0.0, invalid_value);

    //+ If the parameter types do not match, an error message is
    //+ printed and the value 0 of the requested type is returned.
    int invalid = params.get_value<int>("PI");
    EXPECT_EQ(0, invalid);
    
    //+ If the parameter types do not match, an error message is
    //+ printed and an empty vector of the requested type is returned.
    std::vector<double> pies = params.get_value<std::vector<double> >("PI");
    EXPECT_EQ(0u, pies.size());
    //ENDParametersErrors

    //+ The parameter query by name is case-sensitive
    double not_found = params.get_value<double>("pi");
    EXPECT_EQ(0.0, not_found);
    
    //BEGINUnsupportedTypes
    //+ Adding a parameter of "unsupported" type...
    stk::util::ParameterList more_params;
    std::complex<double> phase(3.14,2.718);
    more_params.set_param("phase", phase);

    //+ The print system doesn't know about this type, so will print
    //+ a warning message about unrecognized type.
    more_params.write_parameter_list(std::cout);

    //+ However, you can still retrieve the value of the parameter
    //+ if you know what type it is.
    std::complex<double> my_phase = more_params.get_value<std::complex<double> >("phase");
    EXPECT_EQ(my_phase, phase);

    //+ The Parameter class won't help you on determining the type,
    //+ You must know what it is.
    EXPECT_EQ(more_params.get_param("phase").type, stk::util::ParameterType::INVALID);

    //+ If the wrong type is specified, an exception will be thrown...
    EXPECT_ANY_THROW(more_params.get_value<std::complex<int> >("phase"));
    //ENDUnsupportedTypes
  }
}

#endif //STK_HAVE_BOOST

