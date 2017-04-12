#ifndef STK_UTIL_ENVIRONMENT_COMMANDLINEPARSER_HPP
#define STK_UTIL_ENVIRONMENT_COMMANDLINEPARSER_HPP

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

namespace stk {

class CommandLineParser
{
public:
    CommandLineParser() : optionsDesc("Test command line") {}

    template <typename ValueType>
    void add_option(const std::string &option, const std::string &description)
    {
        optionsDesc.add_options()
          (option.c_str(), boost::program_options::value<ValueType>(), description.c_str());
    }

    template <typename ValueType>
    void add_option_with_default(const std::string &option, const std::string &description, const ValueType &defaultValue)
    {
        optionsDesc.add_options()
          (option.c_str(), boost::program_options::value<ValueType>()->default_value(defaultValue), description.c_str());
    }

    void parse(int argc, const char *argv[])
    {
        try
        {
            boost::program_options::store(boost::program_options::parse_command_line(argc, argv, optionsDesc), varMap);
        }
        catch(std::exception &e)
        {
            std::cerr << "Error: " << e.what() << std::endl;
            throw e;
        }
    }

    bool is_option_provided(const std::string &option)
    {
        return varMap.count(option) > 0;
    }

    bool is_empty()
    {
        return varMap.empty();
    }

    template <typename ValueType>
    ValueType get_option_value(const std::string &option)
    {
        return varMap[option].as<ValueType>();
    }

private:
    boost::program_options::variables_map varMap;
    boost::program_options::options_description optionsDesc;
};

}

#endif //STK_UTIL_ENVIRONMENT_COMMANDLINEPARSER_HPP
