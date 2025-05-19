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

#include "stk_util/environment/ParseCommandLineArgs.hpp"
#include "stk_util/environment/OptionsSpecification.hpp"  // for Option, OptionsSpecification
#include "stk_util/environment/ParsedOptions.hpp"         // for ParsedOptions
#include "stk_util/util/ReportHandler.hpp"                // for ThrowRequireMsg
#include "stk_util/util/string_utils.hpp"                 // for dash_it, rm_dashes
#include <cstddef>                                        // for size_t
#include <algorithm>                                      // for fill
#include <memory>                                         // for __shared_ptr_access, shared_ptr
#include <ostream>                                        // for operator<<, basic_ostream
#include <stdexcept>                                      // for runtime_error
#include <string>                                         // for string, operator<<, char_traits
#include <vector>                                         // for vector, vector<>::reference

namespace stk {

bool insert_option(const Option& option, const std::string& value, stk::ParsedOptions& parsedOptions)
{
    std::string theValue = value.empty() ? option.defaultValue : value;
    option.store(theValue);
//std::cerr<<"inserting name="<<option.name<<", abbrev="<<option.abbrev<<", value="<<theValue<<std::endl;
    return parsedOptions.insert(option.name, option.abbrev, theValue);
}

void make_sure_all_required_were_found(const OptionsSpecification& optionsDesc,
                                       const ParsedOptions& parsedOptions)
{
  for(const auto& option : optionsDesc.get_options_plus_sub_options()) {
    STK_ThrowRequireMsg(!option->isRequired || parsedOptions.count(option->name) == 1,
                    "Error in stk::parse_command_line_args: "
                    << "Required option '"<<dash_it(option->name)<<"' not found");
  }
}

void parse_command_line_args(int argc, const char** argv,
                             const OptionsSpecification& optionsDesc,
                             stk::ParsedOptions& parsedOptions)
{
  for(const auto& option : optionsDesc.get_options_plus_sub_options()) {
    if (!option->defaultValue.empty() && !option->isImplicit) {
      insert_option(*option, "", parsedOptions);
    }
  }

  std::vector<bool> argHasBeenUsed(argc, false);
  std::vector<std::pair<const char*,int>> posArgs;
  posArgs.reserve(argc);
  for(int i=1; i<argc; ++i) {
    const char* arg = argv[i];
    if (arg == nullptr) {
      continue;
    }

    const bool isPositionalArg = (arg[0] != '-');
    if (isPositionalArg) {
      posArgs.push_back(std::make_pair(arg,i));
      continue;
    }

    std::string optionKey = rm_dashes(std::string(arg));
    std::string value;
    bool argContainedEquals = false;

    size_t posOfEqualSign = optionKey.find("=");
    if (posOfEqualSign != std::string::npos) {
      argContainedEquals = true;
      value = optionKey.substr(posOfEqualSign+1);
      optionKey = optionKey.substr(0, posOfEqualSign);
      STK_ThrowRequireMsg(!value.empty(), "stk::parse_command_line_args: Missing value for option " << dash_it(optionKey));
    }

    const Option& option = optionsDesc.find_option(optionKey);
    if (option.name.empty()) {
      continue;
    }

    argHasBeenUsed[i] = true;

    const bool useNextArgAsValue = !option.isFlag && !argContainedEquals;

    if (useNextArgAsValue) {
      const bool nextArgContainsAValue = ((argc >= 3) && (i < (argc-1)) && (argv[i+1] != nullptr) && (argv[i+1][0] != '-'));
      if (nextArgContainsAValue) {
        value = argv[i+1];
        argHasBeenUsed[i+1] = true;
        ++i;
      }
      else {
        STK_ThrowRequireMsg(option.isImplicit, "stk::parse_command_line_args: Missing value for option " << dash_it(option.name));
        value = option.defaultValue;
      }
    }

    insert_option(option, value, parsedOptions);
    parsedOptions.set_parsed(option.name);
  }

  const int numExpectedPositionalOptions = optionsDesc.get_num_positional_options();

  bool foundRedundantFlagAndPositionalArgs = false;
  std::string errMsg;
  if (numExpectedPositionalOptions > 0) {
    int positionalArgIndex = 0;

    for(int posIdx = 0; posIdx < numExpectedPositionalOptions; ++posIdx) {
      const Option& option = optionsDesc.get_positional_option(posIdx);
      const bool alreadyParsed = parsedOptions.is_parsed(option.name);

      if (!alreadyParsed && static_cast<unsigned>(positionalArgIndex) < posArgs.size()) {
        const char* arg = posArgs[positionalArgIndex].first;
        const int argIdx = posArgs[positionalArgIndex].second;

        insert_option(option, std::string(arg), parsedOptions);

        argHasBeenUsed[argIdx] = true;
        ++positionalArgIndex;
      }
    }

    if (static_cast<unsigned>(positionalArgIndex) != posArgs.size()) {
      foundRedundantFlagAndPositionalArgs = true;
      errMsg = std::string("parse_command_line_args: Detected that 1 or more arguments was redundantly specified (with both a flag and a positional argument.");
    }
  }

  if (parsedOptions.count("help") == 0 && parsedOptions.count("version") == 0) {
    make_sure_all_required_were_found(optionsDesc, parsedOptions);
  }

  if (optionsDesc.is_error_on_unrecognized()) {
    for(int i=1; i<argc; ++i) {
      if (!argHasBeenUsed[i]) {
        throw std::runtime_error(std::string("Unrecognized option: '") + std::string(argv[i]) + "'");
      }
    }
  }

  STK_ThrowErrorMsgIf(foundRedundantFlagAndPositionalArgs, errMsg);
}

}

