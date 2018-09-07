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
// 

#ifndef STK_UTIL_SIERRA_TRACE_HPP
#define STK_UTIL_SIERRA_TRACE_HPP

#include <stddef.h>                     // for size_t
#include <algorithm>                    // for copy
#include <cstring>                      // for strcpy
#include <exception>                    // for uncaught_exception
#include <map>                          // for map, map<>::value_compare
#include <ostream>                      // for ostream
#include <stk_util/util/Writer_fwd.hpp> // for LogMask::LOG_TRACE, etc
#include <string>                       // for string
#include <vector>                       // for vector, vector<>::iterator
namespace stk { namespace diag { class Writer; } }

#define SLIB_TRACE_COVERAGE

namespace stk {
namespace diag {

///
/// @addtogroup DiagTraceDetail
/// @{
///

/**
 * @brief Class <b>Tracespec</b> dissects file specification strings.  It
 * contains a single <b>char</b> const pointer to a function specification string.
 * Accessor functions can dissect the function specification string and return various
 * components of it.
 *
 * @ref DiagTracingDetail
 */
class Tracespec
{
public:
  /**
   * Creates a new <b>Tracespec</b> instance.
   *
   * @param function_spec  a <b>char</b> const pointer to the function full
   *        specification.
   */
  explicit Tracespec(const char *function_spec)
    : m_functionSpec(function_spec)
  {}

  /**
   * @brief Member function <b>getFunctionName</b> returns the function's name.
   *
   * @return      a <b>char</b> const pointer to the function's name.
   */
  std::string getFunctionName() const;

protected:
  const char *      m_functionSpec;    ///< The member function specification
};

typedef std::vector<const char *> TracebackStack;

/**
 * @brief Class <b>Traceback</b> is a stack of <b>char</b> constant
 * pointers to function specifications which have been encounter by Traceback during
 * exception stack unwind.
 *
 * @ref DiagTracingDetail
 */
class Traceback : public Tracespec
{
public:
  /**
   * @brief Enumeration <b>TracebackState</b> lists the traceback execution states.
   *
   * <UL>
   * <LI>RUNNING - No throws are in progress.
   * <LI>THROWING - A throw is in progress.
   * </UL>
   */
  enum TracebackState{RUNNING, THROWING};

  enum {STACK_SIZE = 2048};        ///< Running function spec stack size
  typedef const char *Stack[STACK_SIZE] ;    ///< Stack type
  typedef const char * const *  const_stack_iterator;  ///< const iterator thru stack
  typedef const char **  stack_iterator;      ///< iterator thru stack

  /**
   * @brief Typedef <b>Coverage</b> declares the function usage coverage data type.
   *
   * NOTE: This does really want to use const char * as the key.  I know, things could
   * move, but the definition of a function spec is that it is a string literal only.  So,
   * let it go.
   */
  struct StringLiteralLess
  {
    inline bool operator()(const char *lhs, const char *rhs) const
    {
      return std::strcmp(lhs, rhs) < 0;
    }
  };
  typedef std::map<const char *, int, StringLiteralLess> Coverage;

  /**
   * @brief Class <b>Traceback::Preserve</b> serves as a sentry for traceback
   * stack preservation during additional extension and unwinding.
   */
  class Preserve
  {
  public:
    /**
     * @brief Creates a new <b>Traceback::Preserve</b> sentry.  When the sentry
     * is in place, the traceback stack is preserved from extension and clearing.
     */
    Preserve() {
      Traceback::preserveStack();
    }

    /**
     * @brief Destroys a <b>Preserve</b> sentry which allows traceback stack
     */
    ~Preserve() {
      Traceback::releaseStack();
    }
  };

  /**
   * @brief Creates a new <b>Trace</b> instance, resulting in the printing of
   * the member function name and pushing the depth.
   *
   * @param function_spec  a <b>char</b> variable ...
   */
  explicit Traceback(const char *function_spec)
    : Tracespec(function_spec)
  {
    if (s_top >= &s_stack[STACK_SIZE - 1] || s_top == nullptr)
      s_top = s_stack;
    *s_top++ = function_spec;

    if (s_tracebackState == THROWING && !s_tracebackPreserve && !std::uncaught_exception())
      s_tracebackState = RUNNING;

#ifdef SLIB_TRACE_COVERAGE
    if (s_coverageEnabled)
      ++s_coverage[function_spec];
#endif
  }

  /**
   * @brief Destroys a <b>Traceback</b> instance, resulting in the pushing of
   * the function specification if unwinding the stack.
   */
  ~Traceback() {
    if (!s_tracebackPreserve && std::uncaught_exception() && s_tracebackState == RUNNING) {
      s_tracebackState = THROWING;
      s_storedTop = s_storedStack + (s_top - s_stack);
      std::copy(s_stack, s_top, s_storedStack);
    }
    if (s_top > &s_stack[0])
      --s_top;
  }

  static TracebackStack snapshot();

  /**
   * @brief Member function <b>preserveStack</b> increments the traceback stack
   * preservation counter.
   */
  inline static void preserveStack() {
    ++s_tracebackPreserve;
  }

  /**
   * @brief Member function <b>releaseStack</b> decrements the traceback stack
   * preservation counter.
   */
  inline static void releaseStack() {
    --s_tracebackPreserve;
  }

  /**
   * @brief Member function <b>coverageEnabled</b> returns true if coverage has been
   * enabled.
   *
   * @return			a <b>bool</b> value of true if coverage has been
   *				enabled.
   */
  inline static bool coverageEnabled() {
    return s_coverageEnabled;
  }

  /**
   * @brief Member function <b>getTracebackState</b> returns the value of the
   * traceback state.
   */
  inline static TracebackState getTracebackState() {
    return s_tracebackState;
  }

  /**
   * @brief Class <b>PrintCoverage</b> is a type holder class for printing the stack.
   */
  struct PrintCoverage
  {};

  /**
   * @brief Member function <b>printCoverage</b> creates a
   * <b>PrintCoverage</b> type holder class which enables
   * <b>operator&lt;&lt;</b> to put an coverage to an ostream.
   *
   * @return			a <b>PrintCoverage</b> object.
   */
  static PrintCoverage printCoverage() {
    return PrintCoverage();
  }

  /**
   * @brief Member function <b>printCoverage</b> ...
   *
   * @param os		a <b>std::ostream</b> variable ...
   * @return a <b>std::ostream</b> ...
   */
  static std::ostream &printCoverage(std::ostream &os);

  /**
   * @brief Member function <b>printTraceback</b> writes the traceback stack
   * function specifications to the output stream <b>os</b>.
   *
   * @return      a <b>std::ostream</b> reference to the output stream.
   */
  static std::string printTraceback(const TracebackStack &traceback_stack);

private:
  static TracebackState    s_tracebackState;  ///< State of the traceback system
  static int      s_tracebackPreserve;  ///< Preserve traceback stack
  static int      s_tracebackDisplay;  ///< Display traceback stack
  static const char **    s_top;      ///< Pointer to the top + 1 of the stack
  static Stack      s_stack;    ///< Running functionspec stack
  static const char **    s_storedTop;    ///< Pointer to the top + 1 of the stored stack
  static Stack      s_storedStack;    ///< Stored functionspec stored stack
  static bool			s_coverageEnabled;	///< Coverage is enabled
  static Coverage		s_coverage;		///< Function usage coverage
};


/**
 * @brief Class <b>Trace</b> serves as a sentry for entering routines.  Creating a
 * trace object prints the specified member function name to the specified diag_writer and
 * pushes the diag_writer depth.  On destruction, it prints the member function name again
 * and pops the depth.
 *
 * A tracing depth feature has been incorporated into the Writer which enables
 * diagnostic output only when tracing is activated using a specific
 *
 * @ref DiagTracingDetail
 */
class Trace : public Traceback
{
public:
  /**
   * @brief Typedef <b>ExtraFuncPtr</b> declares the extra function pointer
   * signature.
   */
  typedef Writer & (*ExtraFuncPtr)(Writer &);

  /**
   * @brief Typedef <b>TraceList</b> declares the trace list data type.
   */
  struct TraceList : public std::vector<const char *> {
  public:
    TraceList() {
      s_traceListExists = true;
    }

    ~TraceList() {
      for(size_t i=0; i<size(); ++i) {
        delete [] operator[](i);
      }
      s_traceListExists = false;
    }
  };

  /**
   * @brief Enumeration to describe the trace back flags.
   */
  enum {
    IN_TRACE_LIST  = 0x01    ///< Tracing enabled by this function
  };

  /**
   * @brief Creates a new <b>Trace</b> instance, resulting in the printing of
   * the member function name and pushing the depth.
   *
   * @param dout    a <b>Writer</b> reference to the diagnostic
   *        writer to send trace messages to.
   *
   * @param function_name  a <b>char</b> const pointer to the function
   *        specification. THIS POINTER MUST CONTINUE TO EXIST.
   *
   * @param print_mask    an <b>int</b> value of the diagnostic writer print
   *        mask to enable tracing.
   *
   * @param do_trace            a <b>bool</b> that provide an extra dynamic means of
   *                            turning off tracing.
   */
  Trace(Writer &dout, const char *function_name, int print_mask = LOG_TRACE, bool do_trace = true);

  /**
   * @brief Destroys a <b>Trace</b> instance, resulting in the printing of the
   * member function name and popping the diag_writer depth.
   */
  ~Trace();

  /**
   * @brief Member function <b>setExtra</b> sets the extra function which is called
   * during each trace construction and destrution.  (Not implemented)
   *
   * @param extra    an <b>ExtraFuncPtr</b> to the new extra function.
   *
   * @return      an <b>ExtraFuncPtr</b> to the previous extra
   *        function.
   */
  inline static ExtraFuncPtr setExtra(ExtraFuncPtr extra) {
    ExtraFuncPtr x = s_extra;
    s_extra = extra;
    return x;
  }

  /**
   * @brief Member function <b>addTraceFunction</b> adds a function prefix to the
   * list of function prefixes search to enable tracing.
   *
   * @param function_prefix  a <b>std::string</b> const reference to the function
   *        signature prefix.
   */
  inline static void addTraceFunction(const std::string &function_prefix) {
    char *s = std::strcpy(new char[function_prefix.length() + 1], function_prefix.c_str());

    s_traceList.push_back(s);
  }

  /**
   * @brief Member function <b>clearTraceFunctions</b> removes all function prefixes
   * from the function signature prefix list.
   */
  inline static void clearTraceFunctions() {
    for (auto it = s_traceList.begin(); it != s_traceList.end(); ++it) {
      delete[] (*it);
    }

    s_traceList.clear();
  }

  /**
   * @brief Member function <b>dump</b> writes the trace to the specified
   * Writer.
   *
   * @param dout    a <b>Writer</b> variable reference to write the trace to.
   *
   * @return      a <b>Writer</b> reference to <it>Writer</it>.
   */
  Writer &verbose_print(Writer &dout) const;

private:
  Writer &   m_diagWriter;    ///< The diagnostic writer to write to.
  double     m_startCpuTime;  ///< Cpu time at trace start
  PrintMask  m_lineMask;      ///< Line mask of trace
  bool       m_do_trace;      ///< Force trace off
  int        m_flags;         ///< flags if trace enabled

  static ExtraFuncPtr    s_extra;   ///< Extra output on trace
  static TraceList    s_traceList;  ///< List of functions to trace
  static bool   s_traceListExists;  ///< Flag indicating that the tracelist has been constructed and not destructed
};

/**
 * @brief Member function <b>operator&lt;&lt;</b> writes the trace data to the
 * diagnostic writer.
 *
 * @param dout      a <b>Writer</b> reference to the diagnostic
 *        writer to write to.
 *
 * @param diag_trace    a <b>Trace</b> const reference to the trace
 *        object to write.
 *
 * @return      a <b>Writer</b> reference to the diagnostic
 *        writer.
 */
inline Writer &operator<<(Writer &dout, const Trace &diag_trace) {
  return diag_trace.verbose_print(dout);
}

/**
 * @brief Member function <b>operator&lt;&lt;</b> writes the coverage to the output
 * stream.
 *
 * @param os			a <b>std::ostream</b> reference to the output stream
 *				to write to.
 *
 * @return			a <b>std::ostream</b> reference to the output
 *				stream.
 */
inline std::ostream &operator<<(std::ostream &os, const Traceback::PrintCoverage &) {
  return Traceback::printCoverage(os);
}

} // namespace diag
} // namespace stk

namespace sierra {
namespace Diag {

typedef stk::diag::Tracespec Tracespec;
typedef stk::diag::Traceback Traceback;
typedef stk::diag::Trace Trace;

} // namespace Diag 
} // namespace sierra 


///
/// @}
///

#endif // STK_UTIL_SIERRA_TRACE_HPP
