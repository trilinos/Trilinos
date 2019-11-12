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

#include <stk_util/stk_config.h>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/command_line/CommandLineParserUtils.hpp>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <map>
#include <vector>

namespace stk {

struct Option
{
  std::string name;
  std::string abbrev;
  std::string description;
  std::string defaultValue;
  bool isFlag;
  bool isRequired;
  int position;
};

class OptionsSpecification
{
public:
  OptionsSpecification(const std::string& usagePreamble_="")
   : usagePreamble(usagePreamble_), options() {}
  ~OptionsSpecification(){}

  bool empty() const { return options.empty(); }

  size_t size() const { return options.size(); }

  const Option& find_option(const std::string& nameOrAbbrev) const
  {
    static Option emptyOption;
    unsigned numPartialMatches = 0;
    unsigned idxOfPartialMatch = 0;
    unsigned idx = 0;
    for(const Option& option : options) {
      if (option.name == nameOrAbbrev || option.abbrev == nameOrAbbrev) {
        return option;
      }
      else if ((nameOrAbbrev.size() > 1) && (nameOrAbbrev.size() < option.name.size()) &&
               (option.name.find(nameOrAbbrev)==0)) {
        ++numPartialMatches;
        idxOfPartialMatch = idx;
      }
      ++idx;
    }
    //if there were no complete matches, return a partial match if there was exactly one.
    if (numPartialMatches == 1) {
      return options[idxOfPartialMatch];
    }
    return emptyOption;
  }

  size_t get_num_positional_options() const {
    size_t numPositionalOptions = 0;
    for(const Option& option : options) {
      if (option.position >= -1) {
        ++numPositionalOptions;
      }
    }
    return numPositionalOptions;
  }

  const Option& get_positional_option(int position) const
  {
    for(const Option& option : options) {
      if (option.position == position || option.position == -1) {
        return option;
      }
    }
    ThrowRequireMsg(false, "stk::OptionsSpecification failed to find positional option "<<position);
    static Option emptyOption;
    return emptyOption;
  }

  const std::vector<Option>& get_options() const { return options; }

  OptionsSpecification& add_options() { return *this; }

  OptionsSpecification& operator()(const std::string& spec, const std::string& description)
  {
    std::string defaultValue("");
    return add_option(spec, description, defaultValue, true, false, -2);
  }

  OptionsSpecification& operator()(const std::string& spec,
                                 bool isFlag,
                                 bool isRequired,
                                 const std::string& description,
                                 int position=-2)
  {
    std::string defaultValue("");
    return add_option(spec, description, defaultValue, isFlag, isRequired, position);
  }

  template<typename T>
  OptionsSpecification& operator()(const std::string& spec,
                                 const std::string& description,
                                 const T& defaultValue,
                                 bool isFlag = false,
                                 bool isRequired = false,
                                 int position = -2)
  {
    return add_option(spec, description, defaultValue, isFlag, isRequired, position);
  }

  friend std::ostream& operator<<(std::ostream& out, const OptionsSpecification& od);

private:
  template<typename T>
  OptionsSpecification& add_option(const std::string& spec,
                                 const std::string& description,
                                 const T& defaultValue,
                                 bool isFlag,
                                 bool isRequired,
                                 int position)
  {
    std::string optionName = get_substring_before_comma(spec);
    std::string optionAbbrev = get_substring_after_comma(spec);
    std::ostringstream oss;
    oss<<defaultValue;
    std::string defaultValueStr = oss.str();

    for(const Option& opt : options) {
      if (opt.name == optionName) {
        return *this;
      }
    }

    Option option{optionName, optionAbbrev, description, defaultValueStr,
                  isFlag, isRequired, position};
    options.push_back(option);

    return *this;
  }

  std::string usagePreamble;
  std::vector<Option> options;
};

inline
std::ostream& operator<<(std::ostream& out, const OptionsSpecification& od)
{
  if (!od.usagePreamble.empty()) {
    out << od.usagePreamble << std::endl << std::endl;
  }
  int maxOptStrLen = 0;
  for(const Option& opt : od.options) {
    std::string optionString(dash_it(opt.name)+(!opt.abbrev.empty() ? ",-"+opt.abbrev : ""));
    int strLen = optionString.size();
    maxOptStrLen = std::max(maxOptStrLen,strLen);
  }
  int padWidth = maxOptStrLen+3;
  for(const Option& opt : od.options) {
    std::string optionString(dash_it(opt.name)+(!opt.abbrev.empty() ? ",-"+opt.abbrev : ""));
    out << optionString<<std::setw(padWidth+opt.description.size()-optionString.size())<<opt.description<<(opt.isRequired ? " (required)":"")
        <<(!opt.defaultValue.empty() ? (" default: "+opt.defaultValue) : "")
        <<std::endl;
  }

  return out;
}

}

#endif //STK_UTIL_COMMANDLINE_OPTIONSSPECIFICATION_HPP
