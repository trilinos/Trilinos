/*
 * IOParameterList.h
 *
 *  Created on: Oct 2, 2013
 *      Author: swbova
 */

#ifndef PARAMETERLIST_H_
#define PARAMETERLIST_H_

#include <boost/any.hpp>
#include <map>
#include <vector>
#include <string>
#include <iostream>

namespace stk {
  namespace util {

    namespace ParameterType{
      enum Type {
	INVALID = 0,
	FLOAT,
	DOUBLE,
	INTEGER,
	INT64,
	STRING,
	FLOATVECTOR,
	DOUBLEVECTOR,
	INTEGERVECTOR,
	INT64VECTOR,
	STRINGVECTOR
      };

      namespace {
	template <typename T>
	ParameterType::Type get_type(const T &value)
	{ return ParameterType::INVALID; }

	template <> ParameterType::Type get_type(const double &value)
	{ return ParameterType::DOUBLE; }

	template <> ParameterType::Type get_type(const float &value)
	{ return ParameterType::FLOAT; }

	template <> ParameterType::Type get_type(const int &value)
	{ return ParameterType::INTEGER; }

	template <> ParameterType::Type get_type(const int64_t &value)
	{ return ParameterType::INT64; }

	template <> ParameterType::Type get_type(const std::string &value)
	{ return ParameterType::STRING; }

	template <>
	ParameterType::Type get_type(const std::vector<float> &value)
	{ return ParameterType::FLOATVECTOR; }

	template <>
	ParameterType::Type get_type(const std::vector<double> &value)
	{ return ParameterType::DOUBLEVECTOR; }

	template <>
	ParameterType::Type get_type(const std::vector<int> &value)
	{return ParameterType::INTEGERVECTOR; }

	template <>
	ParameterType::Type get_type(const std::vector<int64_t> &value)
	{ return ParameterType::INT64VECTOR; }

	template <>
	ParameterType::Type get_type(const std::vector<std::string> &value)
	{ return ParameterType::STRINGVECTOR; }
      }
    }

    /*
     * can call set_value without a prior call to set_param, then it
     * just uses the defult "non-value" Param member data
     */
    struct Parameter{
      boost::any value;
      ParameterType::Type type;
      bool toResultsFile;
      bool toRestartFile;

      Parameter()
	: value(),
	  type(ParameterType::INVALID),
	  toResultsFile(false),
	  toRestartFile(false)
      {}
    };

    typedef std::map<const std::string, Parameter> ParameterMapType;

    class ParameterList {
    public:
      template <typename T> void set_param (const std::string &name,
					    const T value,
					    bool toOutput = false,
					    bool toRestart = false)
      {
        Parameter p;
        p.value = value;
	p.type = ParameterType::get_type(value);
        p.toResultsFile = toOutput;
        p.toRestartFile = toRestart;
        parameterData[name] = p;
      }

      template <typename T> void set_value (const std::string &name, const T value)
      {
        Parameter &p = parameterData[name];
	p.type = ParameterType::get_type(value);
	if (p.type == ParameterType::INVALID) {
	  std::cerr << "WARNING: Parameter named '" << name
		    << "' is not a supported type.\n";
	}
	p.value = value;
      }

      Parameter get_param (const std::string name)
      {
        ParameterMapType::iterator it = parameterData.find(name);
        if(it == parameterData.end() ) {
	  std::cerr << "ERROR: Parameter named '" << name << "' not found\n";
        }
        return (*it).second;
      }

      template <typename T> T get_value (const std::string name)
      {
        Parameter p = get_param(name);
	if (p.type == ParameterType::get_type(T())) {
	  return boost::any_cast<T>(p.value);
	} else {
	  std::cerr << "ERROR: Parameter named '" << name
		    << "' has an incorrect type specfied for the get_value"
		    << " template type.\n";
	}
	return T();
      }

      void write_parameter_list(std::ostream & stream);
      ParameterMapType::const_iterator begin() { return parameterData.begin();}
      ParameterMapType::const_iterator end() { return parameterData.end();}

    private:
      ParameterMapType parameterData;

    };
  }
}


#endif /* PARAMETERLIST_H_ */
