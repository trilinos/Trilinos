#include <stk_util/util/ParameterList.hpp>

namespace stk {
  namespace util {

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
			  << " and the values are ";
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
			  << " and the values are ";
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
			  << " and the values are ";
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
			  << " and the values are ";
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
			  << " and the value is "
			  << boost::any_cast<std::string>(parameter.value)
			  << std::endl;
		break;
	      }
	    case ParameterType::STRINGVECTOR:
	      {
		std::vector<std::string> vec = boost::any_cast<std::vector<std::string> >(parameter.value);
		stream << "Parameter '" << parameterName << "' is of type vector of strings"
			  << " and the values are ";
		for (size_t j = 0; j < vec.size(); ++j) {
		  if (j>0) stream << ", ";
		  stream << vec[j];
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
