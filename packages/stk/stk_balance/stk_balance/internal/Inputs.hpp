// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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
#ifndef TPLS_SRC_TRILINOS_PACKAGES_STK_STK_BALANCE_STK_BALANCE_INTERNAL_INPUTS_HPP_
#define TPLS_SRC_TRILINOS_PACKAGES_STK_STK_BALANCE_STK_BALANCE_INTERNAL_INPUTS_HPP_

#include <string.h>
#include <string>

#include <sys/stat.h> // move us
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>

namespace stk { namespace balance {

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

}
}



#endif /* TPLS_SRC_TRILINOS_PACKAGES_STK_STK_BALANCE_STK_BALANCE_INTERNAL_INPUTS_HPP_ */
