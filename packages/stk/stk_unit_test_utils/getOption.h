#ifndef UNITTESTUTILS_OPTIONS_PARSING
#define UNITTESTUTILS_OPTIONS_PARSING

#include <string>
#include <sstream>

extern int gl_argc;
extern char** gl_argv;

namespace stk
{
namespace unit_test_util
{

inline std::string get_option(const std::string& option, const std::string defaultString="no")
{
    std::string returnValue = defaultString;
    if ( gl_argv != 0 )
    {
        for (int i=0;i<gl_argc;i++)
        {
            std::string input_argv(gl_argv[i]);
            if ( option == input_argv )
            {
                if ( (i+1) < gl_argc )
                {
                    returnValue = std::string(gl_argv[i+1]);
                }
                break;
            }
        }
    }
    return returnValue;
}

template <typename T>
T get_command_line_option(const std::string &option, const T &defaultValue)
{
    std::string str = get_option(option, std::to_string(defaultValue));
    std::istringstream ss(str);
    T val=0;
    ss >> val;
    return val;
}

}
}

#endif

