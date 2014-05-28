/**   ------------------------------------------------------------
 *    Copyright 2003 - 2011 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <string>
#include <cstring>
#include <sstream>
#include <list>
#include <exception>

#include <stk_util/diag/Trace.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/environment/FormatTime.hpp>
#include <stk_util/diag/Platform.hpp>

namespace stk {
namespace diag {

Trace::ExtraFuncPtr
Trace::s_extra = 0;

Trace::TraceList
Trace::s_traceList;

bool
Trace::s_traceListExists = false;

Traceback::stack_iterator
Traceback::s_top = s_stack;

Traceback::stack_iterator
Traceback::s_storedTop = s_storedStack;

Traceback::Stack
Traceback::s_stack;

Traceback::Stack
Traceback::s_storedStack;

Traceback::TracebackState
Traceback::s_tracebackState = Traceback::RUNNING;

int
Traceback::s_tracebackPreserve = 0;

int
Traceback::s_tracebackDisplay = 0;

bool
Traceback::s_coverageEnabled = false;

Traceback::Coverage
Traceback::s_coverage;

namespace {

bool
prefix_compare(
  const char *    prefix,
  const char *    name)
{
  for (; *prefix != 0 && *prefix == *name; ++prefix, ++name)
    ;
  return *prefix == 0;
}


bool
prefix_find(
  const Trace::TraceList &  trace_list,
  const char *    s)
{
  for (Trace::TraceList::const_iterator it = trace_list.begin(); it != trace_list.end(); ++it)
    if (prefix_compare((*it), s))
      return true;;
  return false;
}


size_t
get_heap_used()
{
  return sierra::Env::get_heap_usage();
}

std::string::const_iterator
find_next_char(
  std::string::const_iterator  p,
  std::string::const_iterator  end,
  char        c)
{
  while (p != end && *p != c)
    p++;
  return p;
}

std::string::const_iterator
find_prev_char(
  std::string::const_iterator  begin,
  std::string::const_iterator  p,
  char        c)
{
  while (p != begin && *p != c)
    p--;
  return p;
}

std::string::const_iterator
find_prev_double_colon(
  std::string::const_iterator  begin,
  std::string::const_iterator  p)
{
  std::string::const_iterator it = p - 1;

  while ((it = find_prev_char(begin, it, ':')) != begin) {
    if (*(it - 1) == ':')
      return it + 1;
    --it;
  }

  return it;
}


inline std::string::const_iterator find_next_open_paren(std::string::const_iterator p, std::string::const_iterator end) {
  return find_next_char(p, end, '(');
}

/**
 * Member function <b>get_function_spec_parts</b> parses <b>spec</b> into the
 * components namespace, class, name and argument list.  The name immediately proceeds the
 * open parenthesis of the argument list; the class immediately proceeds the name
 * separated by double colons; the namespace is everything which immediately proceeds the class.
 *
 * @param spec      a <b>std::string</b> const reference to the function
 *        specification.
 *
 * @param namespace_name  a <b>std::string</b> reference to receive the
 *        namespace part of the function specification.
 *
 * @param class_name    a <b>std::string</b>  reference to receive the
 *        class part of the function specification.
 *
 * @param name      a <b>std::string</b> reference to receive the
 *        name part of the function specification.
 *
 * @param arglist    a <b>std::vector</b> reference of strings to receive
 *        the function arguments. (not implemented)
 */
void
get_function_spec_parts(
  const std::string &    spec,
  std::string &      namespace_name,
  std::string &      class_name,
  std::string &      function_name,
  std::vector<std::string> &  arglist)
{
  namespace_name.erase(namespace_name.begin(), namespace_name.end());
  class_name.erase(class_name.begin(), class_name.end());
  function_name.erase(function_name.begin(), function_name.end());
  arglist.erase(arglist.begin(), arglist.end());

  std::string::const_iterator it_paren = find_next_open_paren(spec.begin(), spec.end());
  std::string::const_iterator it_func_name = find_prev_double_colon(spec.begin(), it_paren);
  function_name = std::string(it_func_name, it_paren);
  if (it_func_name != spec.begin()) {
    it_func_name -= 2;
    std::string::const_iterator it_class_name = find_prev_double_colon(spec.begin(), it_func_name);
    class_name = std::string(it_class_name, it_func_name);
    if (it_class_name != spec.begin()) {
      it_class_name -= 2;
      namespace_name = std::string(spec.begin(), it_class_name);
    }
  }
}

std::string
format_memory(
  int      size)
{
  static const char *suffix[] = {" B", " KB", " MB", " GB"};

  char sign = size < 0 ? '-' : '+';

  size = size > 0 ? size : -size;

  int s = size/10240;

  unsigned int i;
  for (i = 0; i < sizeof(suffix); i++) {
    if (s == 0)
      break;
    size /= 1024;
    s /= 1024;
  }

  std::stringstream strout;

  strout << sign << size << suffix[i];

  return strout.str();
}

} // namespace


std::string
Tracespec::getFunctionNamespace() const
{
  std::string namespace_name;
  std::string class_name;
  std::string function_name;
  std::vector<std::string> arglist;

  get_function_spec_parts(m_functionSpec, namespace_name, class_name, function_name, arglist);

  return namespace_name;
}


std::string
Tracespec::getFunctionClass() const
{
  std::string namespace_name;
  std::string class_name;
  std::string function_name;
  std::vector<std::string> arglist;

  get_function_spec_parts(m_functionSpec, namespace_name, class_name, function_name, arglist);

  return namespace_name + "::" + class_name;
}


std::string
Tracespec::getFunctionShortClass() const
{
  std::string namespace_name;
  std::string class_name;
  std::string function_name;
  std::vector<std::string> arglist;

  get_function_spec_parts(m_functionSpec, namespace_name, class_name, function_name, arglist);

  return class_name;
}


std::string
Tracespec::getFunctionName() const
{
  std::string namespace_name;
  std::string class_name;
  std::string function_name;
  std::vector<std::string> arglist;

  get_function_spec_parts(m_functionSpec, namespace_name, class_name, function_name, arglist);

  return namespace_name + "::" + class_name + "::" + function_name;
}


std::string
Tracespec::getFunctionShortName() const
{
  std::string namespace_name;
  std::string class_name;
  std::string function_name;
  std::vector<std::string> arglist;

  get_function_spec_parts(m_functionSpec, namespace_name, class_name, function_name, arglist);

  return function_name;
}


struct CoverageValueSort
{
  int operator()(const std::pair<const char *, int> &s1, const std::pair<const char *, int> &s2) {
    return std::strcmp(s1.first, s2.first) < 0;
  }
};

Trace::Trace(
  Writer &dout,
  const char *function_name,
  int line_mask,
  bool do_trace)
  : Traceback(function_name),
    m_diagWriter(dout),
    m_startCpuTime(0.0),
    m_startMemAlloc(0),
    m_lineMask(line_mask),
    m_do_trace(do_trace),
    m_flags((dout.isTracing()
             || (dout.shouldTrace(m_lineMask)
                 && (s_traceListExists && (s_traceList.empty() || prefix_find(s_traceList, m_functionSpec))))) ? IN_TRACE_LIST : 0)
{
  if (m_do_trace && (m_flags & IN_TRACE_LIST)) {
    m_diagWriter.incTraceDepth();

    m_diagWriter.m(m_lineMask) << m_functionSpec
                               << (std::uncaught_exception() ? " (throw unwinding) " : "")
                               << push << dendl;

    if (dout.shouldPrint(LOG_TRACE_STATS)) {
      m_startCpuTime = sierra::Env::cpu_now();
      m_startMemAlloc = get_heap_used();
    }
  }
}


Trace::~Trace()
{
  if (m_do_trace && (m_flags & IN_TRACE_LIST)) {
    if (m_diagWriter.shouldPrint(LOG_TRACE_STATS)) {
      m_startCpuTime = sierra::Env::cpu_now() - m_startCpuTime;
      m_startMemAlloc = get_heap_used() - m_startMemAlloc;
    }

    if (m_diagWriter.shouldPrint(LOG_TRACE_STATS)) {
      m_diagWriter.m(m_lineMask) << "[" << stk::formatTime(m_startCpuTime)
                                 << "s, " << format_memory(m_startMemAlloc) << "]" << dendl;
    }

    m_diagWriter.m(m_lineMask) << (std::uncaught_exception() ? " (throw unwinding) " : "")
                               << pop << dendl;

    m_diagWriter.decTraceDepth();
  }
}


Writer &
Trace::verbose_print(
  Writer &		dout) const
{
  dout << "Trace, " << m_functionSpec;
  return dout;
}


TracebackStack
Traceback::snapshot()
{
  TracebackStack traceback_stack;
  traceback_stack.reserve(s_top - s_stack);

  if (Traceback::getTracebackState() == Traceback::RUNNING)
    for (const_stack_iterator it = s_top - 1; it >= s_stack; --it)
      traceback_stack.push_back(*it);
  else
    for (const_stack_iterator it = s_storedTop - 1; it >= s_storedStack; --it)
      traceback_stack.push_back(*it);

  return traceback_stack;
}


std::string
Traceback::printTraceback(
  const TracebackStack &  traceback_stack)
{
  std::ostringstream s;
  if (traceback_stack.empty())
    s << "    traceback not available" << std::endl;
  else {
    for (TracebackStack::const_iterator it = traceback_stack.begin(); it != traceback_stack.end(); ++it)
      s << "    from " << (*it) << std::endl;
  }

  return s.str();
}

std::ostream &
Traceback::printCoverage(
  std::ostream &	os)
{
  std::vector<std::pair<const char *, int> >	sorted_list;
  sorted_list.reserve(s_coverage.size());

  for (Coverage::const_iterator it = s_coverage.begin(); it != s_coverage.end(); ++it)
    sorted_list.push_back(std::pair<const char *, int>((*it).first, (*it).second));

  std::sort(sorted_list.begin(), sorted_list.end(), CoverageValueSort());

  for (std::vector<std::pair<const char *, int> >::const_iterator it = sorted_list.begin(); it != sorted_list.end(); ++it)
    os << "<FUNCTION specification=\"" << (*it).first << "\" count=\"" << (*it).second << "\"/>" << std::endl;

  return os;
}


} // namespace diag
} // namespace stk
