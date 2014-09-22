// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <gtest/gtest.h>
#include <stk_util/util/ParameterList.hpp>

inline void validate_parameters_equal_value(const stk::util::Parameter &parameter,
					    const stk::util::Parameter &gold_parameter)
{
  ASSERT_EQ(parameter.type, gold_parameter.type);
  switch(parameter.type)
    {
    case stk::util::ParameterType::INTEGER:
      {
	ASSERT_EQ(boost::any_cast<int>(parameter.value),
		  boost::any_cast<int>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::INT64:
      {
	ASSERT_EQ(boost::any_cast<int64_t>(parameter.value),
		  boost::any_cast<int64_t>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::DOUBLE:
      {
	ASSERT_EQ(boost::any_cast<double>(parameter.value),
		  boost::any_cast<double>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::FLOAT:
      {
	ASSERT_EQ(boost::any_cast<float>(parameter.value),
		  boost::any_cast<float>(gold_parameter.value));
	break;
      }
    case stk::util::ParameterType::DOUBLEVECTOR:
      {
	std::vector<double> vec = boost::any_cast<std::vector<double> >(parameter.value);
	std::vector<double> gvec = boost::any_cast<std::vector<double> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    case stk::util::ParameterType::FLOATVECTOR:
      {
	std::vector<float> vec = boost::any_cast<std::vector<float> >(parameter.value);
	std::vector<float> gvec = boost::any_cast<std::vector<float> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    case stk::util::ParameterType::INTEGERVECTOR:
      {
	std::vector<int> vec = boost::any_cast<std::vector<int> >(parameter.value);
	std::vector<int> gvec = boost::any_cast<std::vector<int> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    case stk::util::ParameterType::INT64VECTOR:
      {
	std::vector<int64_t> vec = boost::any_cast<std::vector<int64_t> >(parameter.value);
	std::vector<int64_t> gvec = boost::any_cast<std::vector<int64_t> >(gold_parameter.value);
	ASSERT_EQ(vec.size(), gvec.size());
	for (size_t j = 0; j < vec.size(); ++j) {
	  ASSERT_EQ(vec[j], gvec[j]);
	}
	break;
      }
    default:
      ASSERT_EQ(1,0) << "Invalid type found in validate_parameters_equal_value";
    }
}
