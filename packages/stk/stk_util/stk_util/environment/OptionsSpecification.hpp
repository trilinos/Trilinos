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

#ifndef STK_UTIL_COMMANDLINE_OPTIONSSPECIFICATION_HPP
#define STK_UTIL_COMMANDLINE_OPTIONSSPECIFICATION_HPP

#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include "stk_util/util/string_utils.hpp"   // for get_substring_after_comma, get_substring_befo...
#include <cstddef>                          // for size_t
#include <algorithm>                        // for max
#include <memory>                           // for shared_ptr, __shared_ptr_access
#include <sstream>                          // for ostream, ostringstream, istringstream, operat...
#include <string>                           // for string, operator==, operator<<, basic_string
#include <type_traits>                      // for is_same
#include <vector>                           // for vector

namespace stk {

constexpr int INVALID_POSITION = -2;

struct Option
{
  Option()
  : name(), abbrev(), description(), defaultValue(),
    isFlag(false), isRequired(false), isImplicit(false),
    position(INVALID_POSITION)
  {}

  Option(const std::string& nm, const std::string& abrv,
         const std::string& desc, const std::string& deflt,
         bool isFlg, bool isReqd, bool isImplct, int pos)
  : name(nm), abbrev(abrv), description(desc), defaultValue(deflt),
    isFlag(isFlg), isRequired(isReqd), isImplicit(isImplct),
    position(pos)
  {}
  virtual ~Option() {}

  std::string name;
  std::string abbrev;
  std::string description;
  std::string defaultValue;
  bool isFlag;
  bool isRequired;
  bool isImplicit;
  int position;

  virtual void store(const std::string& /* val */) const {}
};

template<typename T>
struct OptionT : public Option
{
  OptionT(const std::string& nm, const std::string& abrv,
         const std::string& desc, const std::string& deflt,
         T* target, int pos=INVALID_POSITION)
  : Option(nm, abrv, desc, deflt, false, false, false, pos),
    targetVariable(target)
  {}

  void store(const std::string& val) const override
  {
    if (!val.empty() && targetVariable != nullptr) {
      std::istringstream iss(val);
      T tVal;
      iss >> tVal;
      STK_ThrowRequireMsg(!iss.fail(), "Error in OptionT, failed to store '"<<val<<"'");
      *targetVariable = tVal;
    }
  }

  mutable T* targetVariable;
};

template<typename T>
struct DefaultValue {
  DefaultValue(const T& val) : value(val) {}
  const T& value;
};

template<typename T>
struct ImplicitValue {
  ImplicitValue(const T& val) : value(val) {}
  const T& value;
};

template<typename T>
struct TargetPointer {
  TargetPointer(T* ptr) : target(ptr) {}
  T* target;
};

template<typename T>
struct ValueType {
  typedef T type;
};

class OptionsSpecification
{
public:
  OptionsSpecification(const std::string& usagePreamble_="", unsigned lineLen=80);

  virtual ~OptionsSpecification(){}

  //copy constructor
  OptionsSpecification(const OptionsSpecification& spec);

  //assignment operator
  OptionsSpecification& operator=(const OptionsSpecification& spec);

  //move constructor
  OptionsSpecification(OptionsSpecification&& spec);

  //move assignment operator
  OptionsSpecification& operator=(OptionsSpecification&& spec);

  bool empty() const { return options.empty(); }

  size_t size() const { return options.size(); }

  void set_error_on_unrecognized() { errorOnUnrecognized = true; }
  bool is_error_on_unrecognized() const { return errorOnUnrecognized; }

  const Option& find_option(const std::string& nameOrAbbrev) const;

  size_t get_num_positional_options() const;

  const Option& get_positional_option(int position) const;

  const std::vector<std::shared_ptr<Option>>& get_options() const { return options; }
  const std::vector<std::shared_ptr<Option>>& get_options_plus_sub_options() const { return optionsPlusSubOptions; }

  OptionsSpecification& add_options() { return *this; }

  OptionsSpecification& add(const OptionsSpecification& spec);

  OptionsSpecification& operator()(const std::string& spec, const std::string& description)
  {
    std::string defaultValue("");
    const bool isFlag = true;
    const bool isRequired = false;
    const bool isImplicit = false;
    const int position = INVALID_POSITION;
    return add_option(spec, description, defaultValue,
                      isFlag, isRequired, isImplicit,
                      position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                   const std::string& description,
                                   ValueType<T> /*valueType*/)
  {
    std::string defaultValue("");
    const bool isFlag = false;
    const bool isRequired = false;
    const bool isImplicit = false;
    const int position = INVALID_POSITION;
    return add_option(spec, description, defaultValue,
                      isFlag, isRequired, isImplicit,
                      position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                   const std::string& description,
                                   ImplicitValue<T> implicitValue)
  {
    const bool isFlag = false;
    const bool isRequired = false;
    const bool isImplicit = true;
    const int position = INVALID_POSITION;
    return add_option(spec, description, implicitValue.value,
                      isFlag, isRequired, isImplicit,
                      position);
  }

  OptionsSpecification& operator()(const std::string& spec,
                                 bool isFlag,
                                 bool isRequired,
                                 const std::string& description,
                                 int position=INVALID_POSITION)
  {
    std::string defaultValue("");
    const bool isImplicit = false;
    return add_option(spec, description, defaultValue,
                      isFlag, isRequired, isImplicit,
                      position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                   const std::string& description,
                                   TargetPointer<T> targetToStoreResult)
  {
    const bool isRequired = true;
    T defaultValue = std::is_same<T,std::string>::value ? "" : 0;
    int position = INVALID_POSITION;
    return add_option(spec, description, defaultValue, targetToStoreResult.target, false, isRequired, position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                   const std::string& description,
                                   TargetPointer<T> targetToStoreResult,
                                   DefaultValue<T> defaultValue,
                                   int position = INVALID_POSITION)
  {
    return add_option(spec, description, defaultValue.value, targetToStoreResult.target, false, false, position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                   const std::string& description,
                                   DefaultValue<T> defaultValue,
                                   TargetPointer<T> targetToStoreResult,
                                   int position = INVALID_POSITION)
  {
    return add_option(spec, description, defaultValue.value, targetToStoreResult.target, false, false, position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                 const std::string& description,
                                 DefaultValue<T> defaultValue,
                                 bool isFlag = false,
                                 bool isRequired = false,
                                 int position = INVALID_POSITION)
  {
    const bool isImplicit = false;
    return add_option(spec, description, defaultValue.value,
                      isFlag, isRequired, isImplicit,
                      position);
  }

  friend std::ostream& operator<<(std::ostream& out, const OptionsSpecification& od);

private:
  template<typename T>
  OptionsSpecification& add_option(const std::string& spec,
                                 const std::string& description,
                                 const T& defaultValue,
                                 bool isFlag,
                                 bool isRequired,
                                 bool isImplicit,
                                 int position)
  {
    std::string optionName = get_substring_before_comma(spec);
    std::string optionAbbrev = get_substring_after_comma(spec);
    std::ostringstream oss;
    oss<<defaultValue;
    std::string defaultValueStr = oss.str();

    for(const auto& opt : options) {
      if (opt->name == optionName) {
        return *this;
      }
    }

    options.push_back(std::shared_ptr<Option>(new Option(optionName, optionAbbrev, description,
                                             defaultValueStr,
                                             isFlag, isRequired, isImplicit,
                                             position)));
    optionsPlusSubOptions.push_back(options.back());

    return *this;
  }

  template<typename T>
  OptionsSpecification& add_option(const std::string& spec,
                                 const std::string& description,
                                 const T& defaultValue,
                                 T* targetToStoreResult,
                                 bool /*isFlag*/,
                                 bool /*isRequired*/,
                                 int position)
  {
    std::string optionName = get_substring_before_comma(spec);
    std::string optionAbbrev = get_substring_after_comma(spec);
    std::ostringstream oss;
    oss<<defaultValue;
    std::string defaultValueStr = oss.str();

    for(const auto& opt : options) {
      if (opt->name == optionName) {
        return *this;
      }
    }

    options.push_back(std::shared_ptr<Option>(new OptionT<T>(optionName, optionAbbrev, description,
                                                          defaultValueStr, targetToStoreResult,
                                                          position)));
    optionsPlusSubOptions.push_back(options.back());

    return *this;
  }

  void print(std::ostream& out) const;

  std::string usagePreamble;
  unsigned lineLength;
  std::vector<std::shared_ptr<Option>> options;
  std::vector<OptionsSpecification> subOptionSpecs;
  std::vector<std::shared_ptr<Option>> optionsPlusSubOptions;
  bool errorOnUnrecognized;
};

inline
std::ostream& operator<<(std::ostream& out, const OptionsSpecification& opspec)
{
  opspec.print(out);

  for(const OptionsSpecification& spec : opspec.subOptionSpecs) {
    spec.print(out);
  }

  return out;
}

}

#endif //STK_UTIL_COMMANDLINE_OPTIONSSPECIFICATION_HPP
