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

#include <stk_util/util/ParameterList.hpp>

#ifdef STK_HAVE_BOOST

#include <stddef.h>                     // for size_t
#include "boost/any.hpp"                // for any_cast

namespace stk {
  namespace util {

    Parameter ParameterList::invalid;
    
    void ParameterList::write_parameter_list(std::ostream & stream)
    {
      ParameterMapType::iterator i = parameterData.begin();
      ParameterMapType::iterator iend = parameterData.end();
      for (; i != iend; ++i)
	{
	  std::string parameterName = (*i).first;
	  Parameter &parameter = (*i).second;
	  switch(parameter.type)
	    {
	    case ParameterType::INTEGER:
	      {
		stream << "Parameter '" << parameterName << "' is of type integer"
		       << " and the value is "
		       << boost::any_cast<int>(parameter.value) << std::endl;
		break;
	      }
	    case ParameterType::INT64:
	      {
		stream << "Parameter '" << parameterName << "' is of type int64"
		       << " and the value is "
		       << boost::any_cast<int64_t>(parameter.value) << std::endl;
		break;
	      }
	    case ParameterType::DOUBLE:
	      {
		stream << "Parameter '" << parameterName << "' is of type double"
		       << " and the value is "
		       << boost::any_cast<double>(parameter.value) << std::endl;
		break;
	      }
	    case ParameterType::FLOAT:
	      {
		stream << "Parameter '" << parameterName << "' is of type float"
		       << " and the value is "
		       << boost::any_cast<float>(parameter.value) << std::endl;
		break;
	      }
	    case ParameterType::DOUBLEVECTOR:
	      {
		std::vector<double> vec = boost::any_cast<std::vector<double> >(parameter.value);
		stream << "Parameter '" << parameterName << "' is of type vector of doubles"
		       << " and the " << vec.size() << " values are ";
		for (size_t j = 0; j < vec.size(); ++j) {
		  if (j>0) stream << ", ";
		  stream << vec[j];
		}
		stream << std::endl;
		break;
	      }
	    case ParameterType::FLOATVECTOR:
	      {
		std::vector<float> vec = boost::any_cast<std::vector<float> >(parameter.value);
		stream << "Parameter '" << parameterName << "' is of type vector of floats"
		       << " and the " << vec.size() << " values are ";
		for (size_t j = 0; j < vec.size(); ++j) {
		  if (j>0) stream << ", ";
		  stream << vec[j];
		}
		stream << std::endl;
		break;
	      }
	    case ParameterType::INTEGERVECTOR:
	      {
		std::vector<int> vec = boost::any_cast<std::vector<int> >(parameter.value);
		stream << "Parameter '" << parameterName << "' is of type vector of integers"
		       << " and the " << vec.size() << " values are ";
		for (size_t j = 0; j < vec.size(); ++j) {
		  if (j>0) stream << ", ";
		  stream << vec[j];
		}
		stream << std::endl;
		break;
	      }
	    case ParameterType::INT64VECTOR:
	      {
		std::vector<int64_t> vec = boost::any_cast<std::vector<int64_t> >(parameter.value);
		stream << "Parameter '" << parameterName << "' is of type vector of int64s"
		       << " and the " << vec.size() << " values are ";
		for (size_t j = 0; j < vec.size(); ++j) {
		  if (j>0) stream << ", ";
		  stream << vec[j];
		}
		stream << std::endl;
		break;
	      }
	    case ParameterType::STRING:
	      {
		stream << "Parameter '" << parameterName << "' is of type string"
		       << " and the value is '"
		       << boost::any_cast<std::string>(parameter.value)
		       << "'" << std::endl;
		break;
	      }
	    case ParameterType::STRINGVECTOR:
	      {
		std::vector<std::string> vec = boost::any_cast<std::vector<std::string> >(parameter.value);
		stream << "Parameter '" << parameterName << "' is of type vector of strings"
		       << " and the " << vec.size() << " values are ";
		for (size_t j = 0; j < vec.size(); ++j) {
		  if (j>0) stream << ", ";
		  stream << "'" << vec[j] << "'";
		}
		stream << std::endl;
		break;
	      }
	    default:
	      {
		stream << "WARNING: '" << parameterName
		       << "' is not a supported type. It's value cannot be output."
		       << std::endl;
		break;
	      }
	    }
	}
    }
  }
}

#endif //STK_HAVE_BOOST

