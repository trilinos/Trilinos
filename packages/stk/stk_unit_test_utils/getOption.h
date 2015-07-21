#ifndef UNITTESTUTILS_OPTIONS_PARSING
#define UNITTESTUTILS_OPTIONS_PARSING

#include <string>

extern int gl_argc;
extern char** gl_argv;

namespace unitTestUtils
{

inline std::string getOption(const std::string& option, const std::string defaultString="no")
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

} // end namespace

#endif

