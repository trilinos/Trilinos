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
        optionsDesc.add_options()
          (get_option_spec(option).c_str(), boost::program_options::value<ValueType>()->default_value(defaultValue), option.description.c_str());
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
            boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(optionsDesc).positional(positionalDesc).run(), varMap);
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
