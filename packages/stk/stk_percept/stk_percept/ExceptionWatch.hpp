#ifndef stk_percept_ExceptionWatch_hpp
#define stk_percept_ExceptionWatch_hpp

#include <stdio.h>
#include <exception>

/** This little jewel is a poor man's stack trace.  To use, insert EXCEPTWATCH at the beginning of each function you want to trace:
 *
 * void func() {
 *    EXCEPTWATCH
 *    //...
 *    }
 * see: http://www.velocityreviews.com/forums/t284190-exception-stack-trace.html
 */

namespace stk { 
  namespace percept { 

    class ExceptionWatch {
      int line_;
      char const* pfname_;
    public:
      ExceptionWatch(int line, char const* pfname) : line_(line),
                                                     pfname_(pfname) {}
      ~ExceptionWatch() {
        if(std::uncaught_exception()) {
          //on purpose
          printf("STACKTRACE::ExceptionWatch: line:\t%d\tfile name:\t%s\n", line_, pfname_);
        }
      }
    };
  }
}

#define TOKENPASTE_LOC(x,y) x ## y
#define TOKENPASTE2_LOC(x,y) TOKENPASTE_LOC(x,y)

#ifndef NDEBUG
#define EXCEPTWATCH stk::percept::ExceptionWatch TOKENPASTE2_LOC(exception_watch_, __COUNTER__ ) (__LINE__, __FILE__)
#else
#define EXCEPTWATCH ((void)0)
#endif

#endif
