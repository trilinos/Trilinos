#ifndef TPLS_SRC_TRILINOS_PACKAGES_STK_STK_BALANCE_STK_BALANCE_INTERNAL_INPUTS_HPP_
#define TPLS_SRC_TRILINOS_PACKAGES_STK_STK_BALANCE_STK_BALANCE_INTERNAL_INPUTS_HPP_

#include <string.h>
#include <string>

#include <sys/stat.h> // move us
#include <cerrno>
#include <cstring>
#include <fstream>

namespace stk { namespace balance {

class Inputs
{
public:
    Inputs(const int argc, const char* inputs[])
    : executableName(inputs[0]), exodusFilename(""), outputDirectory(".")
    {
        if(argc>1)
        {
            exodusFilename = std::string(inputs[1]);
        }

        if(argc>2)
            outputDirectory = std::string(inputs[2]);
    }

    std::string get_executable_name() const { return executableName; }
    std::string get_exodus_filename() const { return exodusFilename; }
    std::string get_output_directory() const { return outputDirectory; }


private:
    std::string executableName;
    std::string exodusFilename;
    std::string outputDirectory;
};

// Only call for proc 0 (or any single proc)
// Copied and modified from Ioss_DatabaseIO.C::create_path
inline
bool create_path(const std::string &path)
{
    bool error_found = false;
    std::ostringstream errmsg;

    const int mode = 0777; // Users umask will be applied to this.

    auto iter = path.begin();
    while(iter != path.end() && !error_found)
    {
        iter = std::find(iter, path.end(), '/');
        std::string path_root = std::string(path.begin(), iter);

        if(iter != path.end())
        {
            ++iter; // Skip past the '/'
        }

        if(path_root.empty())
        { // Path started with '/'
            continue;
        }

        struct stat st;
        if(stat(path_root.c_str(), &st) != 0)
        {
            if(mkdir(path_root.c_str(), mode) != 0 && errno != EEXIST)
            {
                errmsg << "ERROR: Cannot create directory '" << path_root << "' : " << std::strerror(errno) << "\n";
                error_found = true;
            }
        }
        else if(!S_ISDIR(st.st_mode))
        {
            errno = ENOTDIR;
            errmsg << "ERROR: Path '" << path_root << "' is not a directory.\n";
            error_found = true;
        }
    }

    if(error_found)
        std::cerr << errmsg.str();

    return error_found == false;
}

inline
bool should_write_usage_info(const std::string& filename)
{
    return filename.empty();
}

inline
bool does_file_exist(const std::string& filename)
{
    bool exists = true;
    if(!std::ifstream(filename))
        exists = false;
    return exists;
}

}
}



#endif /* TPLS_SRC_TRILINOS_PACKAGES_STK_STK_BALANCE_STK_BALANCE_INTERNAL_INPUTS_HPP_ */
