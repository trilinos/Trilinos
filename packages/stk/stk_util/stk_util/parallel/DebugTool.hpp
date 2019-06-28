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
//

#ifndef stk_util_DebugTool_hpp
#define stk_util_DebugTool_hpp

#include <execinfo.h>
#include <cxxabi.h>
#include <stk_util/stk_config.h>
#include <string>
#include <iostream>

inline std::string demangleFunctionNames(char** symbollist, int addrlen)
{
    std::string mangledNamesString("");
#if defined(__GNUC__) && !defined(__ICC)
    size_t funcnamesize = 256;
    char* funcname = (char*) malloc(funcnamesize);

    for(int i = 1; i < addrlen; i++)
    {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        for(char *p = symbollist[i]; *p; ++p)
        {
            if(*p == '(')
                begin_name = p;
            else if(*p == '+')
                begin_offset = p;
            else if(*p == ')' && begin_offset)
            {
                end_offset = p;
                break;
            }
        }

        if(begin_name && begin_offset && end_offset && begin_name < begin_offset)
        {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            int status;
            char* ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
            if(status == 0)
            {
                mangledNamesString += std::string(ret) + " " + begin_offset + "\n";
            }
            else
            {
                mangledNamesString += std::string(begin_name) + " " + begin_offset + "\n";
            }
        }
        else
        {
            mangledNamesString += symbollist[i] + std::string("\n");
        }
    }

    free(funcname);
#endif
    return mangledNamesString;
}

inline std::string getStackTrace()
{
#if defined(__GNUC__) && !defined(__ICC)
    const int N = 16;
    void *trace[N];
    
    int trace_size = backtrace(trace, N);
    char** mangledFunctionNames = backtrace_symbols(trace, trace_size);

    std::string demangledNames = demangleFunctionNames(mangledFunctionNames, trace_size);
    free(mangledFunctionNames);

    return demangledNames;
#else
    return "none";
#endif
}

#endif

