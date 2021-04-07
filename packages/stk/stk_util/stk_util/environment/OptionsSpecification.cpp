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

#include "stk_util/environment/OptionsSpecification.hpp"
#include "stk_util/util/string_utils.hpp"  // for dash_it, make_vector_of_strings
#include <algorithm>                       // for copy, max
#include <iomanip>                         // for operator<<, setw
#include <utility>                         // for move

namespace stk {

OptionsSpecification::OptionsSpecification(const std::string& usagePreamble_, unsigned lineLen)
 : usagePreamble(usagePreamble_), lineLength(lineLen), options(),
   subOptionSpecs(), optionsPlusSubOptions(),
   errorOnUnrecognized(false)
{
}

//copy constructor
OptionsSpecification::OptionsSpecification(const OptionsSpecification& spec)
 : usagePreamble(spec.usagePreamble),
   lineLength(spec.lineLength),
   options(spec.options),
   subOptionSpecs(spec.subOptionSpecs),
   optionsPlusSubOptions(spec.optionsPlusSubOptions),
   errorOnUnrecognized(spec.errorOnUnrecognized)
{
}

//assignment operator
OptionsSpecification& OptionsSpecification::operator=(const OptionsSpecification& spec)
{
  usagePreamble = spec.usagePreamble;
  lineLength = spec.lineLength;
  options = spec.options;
  subOptionSpecs = spec.subOptionSpecs;
  optionsPlusSubOptions = spec.optionsPlusSubOptions;
  errorOnUnrecognized = spec.errorOnUnrecognized;
  return *this;
}

//move constructor
OptionsSpecification::OptionsSpecification(OptionsSpecification&& spec)
 : usagePreamble(std::move(spec.usagePreamble)),
   lineLength(spec.lineLength),
   options(std::move(spec.options)),
   subOptionSpecs(std::move(spec.subOptionSpecs)),
   optionsPlusSubOptions(std::move(spec.optionsPlusSubOptions)),
   errorOnUnrecognized(spec.errorOnUnrecognized)
{
}

//move assignment operator
OptionsSpecification& OptionsSpecification::operator=(OptionsSpecification&& spec)
{
   usagePreamble = std::move(spec.usagePreamble);
   lineLength = spec.lineLength;
   options = std::move(spec.options);
   subOptionSpecs = std::move(spec.subOptionSpecs);
   optionsPlusSubOptions = std::move(spec.optionsPlusSubOptions);
   errorOnUnrecognized = spec.errorOnUnrecognized;
  return *this;
}

const Option& OptionsSpecification::find_option(const std::string& nameOrAbbrev) const
{
  static Option emptyOption;
  unsigned numPartialMatches = 0;
  unsigned idxOfPartialMatch = 0;
  unsigned idx = 0;
  for(const auto& option : options) {
    if (option->name == nameOrAbbrev || option->abbrev == nameOrAbbrev) {
      return *option;
    }
    else if ((nameOrAbbrev.size() > 1) && (nameOrAbbrev.size() < option->name.size()) &&
             (option->name.find(nameOrAbbrev)==0)) {
      ++numPartialMatches;
      idxOfPartialMatch = idx;
    }
    ++idx;
  }
  //if there were no complete matches, return a partial match if there was exactly one.
  if (numPartialMatches == 1) {
    return *options[idxOfPartialMatch];
  }

  for(const OptionsSpecification& spec : subOptionSpecs) {
    const Option& option = spec.find_option(nameOrAbbrev);
    if (!option.name.empty()) {
      return option;
    }
  }

  return emptyOption;
}

size_t OptionsSpecification::get_num_positional_options() const {
  size_t numPositionalOptions = 0;
  for(const auto& option : options) {
    if (option->position != INVALID_POSITION) {
      ++numPositionalOptions;
    }
  }

  for(const OptionsSpecification& spec : subOptionSpecs) {
    numPositionalOptions += spec.get_num_positional_options();
  }

  return numPositionalOptions;
}

const Option& OptionsSpecification::get_positional_option(int position) const
{
  for(const auto& option : options) {
    if (option->position == position) {
      return *option;
    }
  }

  for(const OptionsSpecification& spec : subOptionSpecs) {
    const Option& option = spec.get_positional_option(position);
    if (!option.name.empty()) {
      return option;
    }
  }

  static Option emptyOption;
  return emptyOption;
}

OptionsSpecification& OptionsSpecification::add(const OptionsSpecification& spec)
{
  subOptionSpecs.push_back(spec);
  optionsPlusSubOptions.insert(optionsPlusSubOptions.end(), spec.options.begin(), spec.options.end());
  return *this;
}

void OptionsSpecification::print(std::ostream& out) const
{
  if (!usagePreamble.empty()) {
    out << usagePreamble << std::endl << std::endl;
  }
  int maxOptStrLen = 0;
  for(const auto& opt : options) {
    std::string optionString(dash_it(opt->name)+(!opt->abbrev.empty() ? ",-"+opt->abbrev : ""));
    int strLen = optionString.size();
    maxOptStrLen = std::max(maxOptStrLen,strLen);
  }

  int padWidth = maxOptStrLen+3;
  int descrWidth = lineLength - padWidth;
  descrWidth = std::max(descrWidth, 20);

  for(const auto& opt : options) {
    std::string optionString(dash_it(opt->name)+(!opt->abbrev.empty() ? ",-"+opt->abbrev : ""));
    std::vector<std::string> descrStrs = make_vector_of_strings(opt->description, ' ', descrWidth);
    out << optionString<<std::setw(padWidth-optionString.size())<<" "<<descrStrs[0]<<std::endl;
    for(unsigned i=1; i<descrStrs.size(); ++i) {
      out << std::setw(padWidth)<<" "<<descrStrs[i]<<std::endl;
    }
    std::string requiredOrDefault = (opt->isRequired ? " (required)":"");
    requiredOrDefault += (!opt->defaultValue.empty() ? (" (default: "+opt->defaultValue+")") : "");
    if (!requiredOrDefault.empty()) {
      out << std::setw(padWidth)<<" "<< requiredOrDefault <<std::endl;
    }
  }

  out << std::endl;
}

}
//namespace stk

