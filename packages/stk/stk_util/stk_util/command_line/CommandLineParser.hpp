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

#ifndef STK_UTIL_ENVIRONMENT_COMMANDLINEPARSER_HPP
#define STK_UTIL_ENVIRONMENT_COMMANDLINEPARSER_HPP

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

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
    enum ParseState { ParseComplete, ParseError, ParseHelpOnly, ParseVersionOnly };
    CommandLineParser() : CommandLineParser("Options") {}
    explicit CommandLineParser(const std::string &usagePreamble) : optionsDesc(usagePreamble)
    {
        add_flag("help,h", "display this help message and exit");
        add_flag("version,v", "display version information and exit");
    }

    void add_flag(const std::string &option, const std::string &description)
    {
        optionsDesc.add_options()
          (option.c_str(), description.c_str());
    }

    template <typename ValueType>
    void add_required_positional(const CommandLineOption &option)
    {
        add_required<ValueType>(option);
        positionalDesc.add(option.name.c_str(), 1);
    }

    template <typename ValueType>
    void add_optional_positional(const CommandLineOption &option, const ValueType &def)
    {
        add_optional<ValueType>(option, def);
        positionalDesc.add(option.name.c_str(), 1);
    }

    template <typename ValueType>
    void add_required(const CommandLineOption &option)
    {
        optionsDesc.add_options()
          (get_option_spec(option).c_str(), boost::program_options::value<ValueType>()->required(), option.description.c_str());
    }

    template <typename ValueType>
    void add_optional(const CommandLineOption &option, const ValueType &defaultValue)
    {
        add_optional(get_option_spec(option), option.description, defaultValue);
    }

    template <typename ValueType>
    void add_optional(const std::string &option, const std::string &description, const ValueType &defaultValue)
    {
        optionsDesc.add_options()
          (option.c_str(), boost::program_options::value<ValueType>()->default_value(defaultValue), description.c_str());
    }

    std::string get_usage() const
    {
        std::ostringstream os;
        os << optionsDesc << std::endl;
        return os.str();
    }

    ParseState parse(int argc, const char *argv[])
    {
        ParseState state = ParseError;
        try
        {
            char** nonconst_argv = const_cast<char**>(argv);
            boost::program_options::store(boost::program_options::command_line_parser(argc, nonconst_argv).options(optionsDesc).positional(positionalDesc).run(), varMap);
            if(is_option_provided("help"))
                return ParseHelpOnly;
            if(is_option_provided("version"))
                return ParseVersionOnly;

            boost::program_options::notify(varMap);
            state = ParseComplete;
        }
        catch(std::exception &e)
        {
            print_message(e.what());
        }
        return state;
    }

    bool is_option_provided(const std::string &option) const
    {
        return varMap.count(option) > 0;
    }

    bool is_empty() const
    {
        return varMap.empty();
    }

    template <typename ValueType>
    ValueType get_option_value(const std::string &option) const
    {
        return varMap[option].as<ValueType>();
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

    boost::program_options::variables_map varMap;
    boost::program_options::options_description optionsDesc;
    boost::program_options::positional_options_description positionalDesc;
};

}

#endif //STK_UTIL_ENVIRONMENT_COMMANDLINEPARSER_HPP
