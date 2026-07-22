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

#ifndef STK_UTIL_COMMAND_LINE_COMMANDLINEPARSER_HPP
#define STK_UTIL_COMMAND_LINE_COMMANDLINEPARSER_HPP

#include "stk_util/environment/OptionsSpecification.hpp"  // for OptionsSpecification, operator<<
#include "stk_util/environment/ParsedOptions.hpp"         // for ParsedOptions
#include "stk_util/util/ReportHandler.hpp"                // for ThrowRequireMsg
#include <iostream>                                       // for endl, operator<<, ostream, ostr...
#include <string>                                         // for string, allocator, operator+

namespace stk {

struct CommandLineOption
{
    std::string name;
    std::string abbreviation;
    std::string description;
};

class CommandLineParser
{
public:
 virtual ~CommandLineParser() = default;
 enum ParseState { ParseComplete, ParseError, ParseHelpOnly, ParseVersionOnly };
 CommandLineParser() : CommandLineParser("Options") {}
 explicit CommandLineParser(const std::string &usagePreamble)
     : optionsSpec(usagePreamble), parsedOptions(), positionalIndex(0)
 {
   add_flag("help,h", "display this help message and exit");
   add_flag("version,v", "display version information and exit");
    }

    void disallow_unrecognized()
    {
      optionsSpec.set_error_on_unrecognized();
    }

    void add_flag(const CommandLineOption &option)
    {
        add_flag(get_option_spec(option), option.description);
    }

    void add_flag(const std::string &option, const std::string &description)
    {
        optionsSpec.add_options()
          (option, description);
    }


    template <typename ValueType>
    void add_required_positional(const CommandLineOption &option)
    {
        add_required_positional<ValueType>(option, positionalIndex);
        ++positionalIndex;
    }

    template <typename ValueType>
    void add_required_positional(const CommandLineOption &option, int position)
    {
        const bool isFlag = false;
        const bool isRequired = true;
        optionsSpec.add_options()
          (get_option_spec(option), isFlag, isRequired, option.description, position);
    }


    template <typename ValueType>
    void add_optional_positional(const CommandLineOption &option, const ValueType &def)
    {
        add_optional_positional<ValueType>(option, def, positionalIndex);
        ++positionalIndex;
    }

    template <typename ValueType>
    void add_optional_positional(const CommandLineOption &option, const ValueType &defaultValue, int position)
    {
        const bool isFlag = false;
        const bool isRequired = false;
        optionsSpec.add_options()
          (get_option_spec(option), option.description, stk::DefaultValue<ValueType>(defaultValue),
           isFlag, isRequired, position);
    }


    template <typename ValueType>
    void add_required(const CommandLineOption &option)
    {
        const bool isFlag = false;
        const bool isRequired = true;
        const int defaultPosition = -2;
        optionsSpec.add_options()
          (get_option_spec(option), isFlag, isRequired, option.description, defaultPosition);
    }


    template <typename ArgValueType>
    void add_optional(const CommandLineOption &option)
    {
        add_optional<ArgValueType>(get_option_spec(option), option.description);
    }

    template <typename ArgValueType>
    void add_optional(const std::string &option, const std::string &description)
    {
        optionsSpec.add_options()
          (option, description, stk::ValueType<ArgValueType>());
    }


    template <typename ValueType>
    void add_optional(const CommandLineOption &option, const ValueType &defaultValue)
    {
        add_optional(get_option_spec(option), option.description, defaultValue);
    }

    template <typename ValueType>
    void add_optional(const std::string &option, const std::string &description, const ValueType &defaultValue)
    {
        const bool isFlag = false;
        const bool isRequired = false;
        const int defaultPosition = -2;
        optionsSpec.add_options()
          (option, description, stk::DefaultValue<ValueType>(defaultValue), isFlag, isRequired, defaultPosition);
    }


    template <typename ValueType>
    void add_optional_implicit(const CommandLineOption &option, const ValueType &defaultValue)
    {
        add_optional_implicit(get_option_spec(option), option.description, defaultValue);
    }

    template <typename ValueType>
    void add_optional_implicit(const std::string &option, const std::string &description, const ValueType &defaultValue)
    {
        optionsSpec.add_options()
          (option, description, stk::ImplicitValue<ValueType>(defaultValue));
    }


    std::string get_usage() const
    {
        std::ostringstream os;
        os << optionsSpec << std::endl;
        return os.str();
    }

    ParseState parse(int argc, const char ** argv);

    bool is_option_provided(const std::string &option) const
    {
        return parsedOptions.count(option) > 0;
    }

    bool is_option_parsed(const std::string& option) const
    {
        return parsedOptions.is_parsed(option);
    }

    bool is_empty() const
    {
        return parsedOptions.empty();
    }

    template <typename ValueType>
    ValueType get_option_value(const std::string &option) const
    {
        STK_ThrowRequireMsg(is_option_provided(option), "Error, option '"<<option<<"'not provided.");
        return parsedOptions[option].as<ValueType>();
    }

protected:
    std::string get_option_spec(const CommandLineOption &option)
    {
        return option.name + "," + option.abbreviation;
    }

    virtual void print_message(const std::string &msg)
    {
        std::cerr << msg << std::endl;
    }

    stk::OptionsSpecification optionsSpec;
    stk::ParsedOptions parsedOptions;
    int positionalIndex;
};

}

#endif //STK_UTIL_COMMAND_LINE_COMMANDLINEPARSER_HPP
