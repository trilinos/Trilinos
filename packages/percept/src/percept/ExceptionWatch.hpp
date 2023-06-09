// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ExceptionWatch_hpp
#define percept_ExceptionWatch_hpp

#include <stdio.h>
#include <exception>

#include <percept/Stacktrace.hpp>


/** This little jewel is a poor man's stack trace.  To use, insert EXCEPTWATCH at the beginning of each function you want to trace:
 *
 * void func() {
 *    EXCEPTWATCH
 *    //...
 *    }
 * see: http://www.velocityreviews.com/forums/t284190-exception-stack-trace.html
 */

  namespace percept {

    class ExceptionWatch {
      int line_;
      char const* pfname_;
    public:
      ExceptionWatch(int line, char const* pfname) : line_(line),
                                                     pfname_(pfname) {}
      ~ExceptionWatch() {
#ifdef __cpp_lib_uncaught_exceptions
        if(std::uncaught_exceptions() > 0) {
#else
        if(std::uncaught_exception()) {
#endif
          //on purpose
          printf("STACKTRACE::ExceptionWatch: line:\t%d\tfile name:\t%s\n", line_, pfname_);
          printf("%s\n", Stacktrace::demangled_stacktrace(30).c_str());
        }
      }
    };
  }

#define TOKENPASTE_LOC(x,y) x ## y
#define TOKENPASTE2_LOC(x,y) TOKENPASTE_LOC(x,y)

#ifndef NDEBUG
#define EXCEPTWATCH percept::ExceptionWatch TOKENPASTE2_LOC(exception_watch_, __COUNTER__ ) (__LINE__, __FILE__)
#else
#define EXCEPTWATCH ((void)0)
#endif

#endif
