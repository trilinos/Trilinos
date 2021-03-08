#pragma once
#ifndef ROL2_TESTFOREXCEPTION_HPP
#define ROL2_TESTFOREXCEPTION_HPP

/** \file  ROL_TestForException.hpp
 *  \brief Defines a macro for conditionally throwing an exception
 *         If the preprocessor variable ROL2_ENABLE_BACKWARD_CPP is
 *         set to true, a stack trace will be performed using backward-cpp
 * 
 *         https://github.com/bombela/backward-cpp
 */

#if ROL2_ENABLE_BACKWARD_CPP // Perform stack trace using 

#include "backward.hpp"

#define ROL_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
{ \
  const bool throw_exception = (throw_exception_test); \
  if(throw_exception) { \
    std::ostringstream omsg; \
    omsg \
      << __FILE__ << ":" << __LINE__ << ":\n\n" \
      << "\n\n" \
      << "Throw test that evaluated to true: "#throw_exception_test \
      << "\n\n" \
      << msg; \
    using namespace backward;  \
    StackTrace st; st.load_here(32); \
    Printer p;  \
    p.object = true; \
    p.color_mode = ColorMode::always; \
    p.address = true;  \
    p.print(st, omsg); \
    const std::string &omsgstr = omsg.str(); \
    throw Exception(omsgstr); \
   } \
}

#else

#define ROL_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
{ \
  const bool throw_exception = (throw_exception_test); \
  if(throw_exception) { \
    std::ostringstream omsg; \
    omsg \
     << __FILE__ << ":" << __LINE__ << ":\n\n" \
     << "\n\n" \
     << "Throw test that evaluated to true: "#throw_exception_test \
     << "\n\n" \
     << msg; \
    const std::string &omsgstr = omsg.str(); \
    throw Exception(omsgstr); \
  }\
}

#endif // ROL2_ENABLE_BACKWARD_CPP








#endif //ROL2_TESTFOREXCEPTION_HPP

