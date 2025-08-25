// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Stacktrace_hpp
#define Stacktrace_hpp

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <cxxabi.h>
#endif

#include <string>
#include <iostream>
#include <sstream>

namespace percept {

#define DO_STACKTRACE_POS 0
// for more efficiency, less info, use the following:
#define DO_STACKTRACE_POS_USE_INT 1

#if DO_STACKTRACE_POS
#  if DO_STACKTRACE_POS_USE_INT
#    define STACKTRACE_POS_I(i) do { Stacktrace::s_position = (i) + __LINE__; } while(0)
#    define STACKTRACE_POS() STACKTRACE_POS_I(0)
#  else
#    define STACKTRACE_POS() do { std::ostringstream str; str << __FILE__ << ":" << __LINE__; Stacktrace::s_position = str.str(); } while(0)
#  endif

#else
#  define STACKTRACE_POS() do { } while(0)
#endif

class Stacktrace {
public:

#if DO_STACKTRACE_POS
#if DO_STACKTRACE_POS_USE_INT
  static int s_position;
#else
  static std::string s_position;
#endif
#endif
  static inline void print_stacktrace(size_t sz=10, const std::string& msg="")
  {
    void **array = new void*[sz];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, sz);

    // print out all the frames to stderr
    fprintf(stderr, "stacktrace: message= %s:\n", msg.c_str());
    //if (!get_rank())
    backtrace_symbols_fd(array, size, 2);
    delete [] array;
  }

  static inline std::string stacktrace(size_t sz=10, const std::string& msg="")
  {
    void **array = new void*[sz];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, sz);

    std::string ret = "stacktrace: "+msg + "\n";
    char ** syms = backtrace_symbols(array, size);
    std::cout << "stacktrace size= " << size << std::endl;
    for (size_t i =0; i < size; ++i)
      {
        ret += std::string(syms[i])+"\n";
      }
    free(syms);
    delete [] array;
    return ret;
  }

  static inline std::string demangle(const std::string& x)
  {

    std::string ret = x;
#if defined(__GNUC__)
    std::string y = x;
    size_t lpos = y.find("(");
    size_t rpos = y.find("+");
    if (lpos != std::string::npos && rpos != std::string::npos)
      {
        y = y.substr(lpos+1,rpos-lpos-1);
      }
    //std::cout << "x= " << x << " x.c_str= " << x.c_str() << " \n y= " << y << std::endl;
    int status=0;
    char *realname=0;
    realname = abi::__cxa_demangle(y.c_str(), 0, 0, &status);
    if (status != 0)
      ret = y;
    else
      {
        ret = realname;
        free(realname);
      }
#endif
    return ret;
  }

  static inline std::string demangled_stacktrace(size_t sz=10, bool also_mangled=false, const std::string& /*msg*/="")
  {
    void **array = new void*[sz];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, sz);

    std::ostringstream str;
#if DO_STACKTRACE_POS
    str << "stacktrace: pos= " << s_position << " " << msg << "\n";
#endif
    std::string ret = str.str();
    char ** syms = backtrace_symbols(array, size);
    std::cout << "stacktrace size= " << size << std::endl;
    if (also_mangled)
      {
        ret += "\n";
        for (size_t i =0; i < size; ++i)
          {
            std::string ss = syms[i];
            ret += ss+"\n";
          }
      }
    ret += "\n";
    for (size_t i =0; i < size; ++i)
      {
        std::string ss = syms[i];
        ret += demangle(ss)+"\n";
      }
    ret += "\n";
    free(syms);
    delete [] array;
    return ret;
  }

};

}//namespace percept

#endif
