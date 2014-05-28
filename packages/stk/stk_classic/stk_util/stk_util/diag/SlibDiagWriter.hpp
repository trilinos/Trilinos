#ifndef STK_UTIL_DIAG_SlibDiagWriter_h
#define STK_UTIL_DIAG_SlibDiagWriter_h

#include <stk_util/diag/Trace.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterParser.hpp>

#include <stk_util/diag/Writer_fwd.hpp>

namespace sierra {
namespace Slib {

stk::diag::Writer &theDiagWriter();

/// Macro <code>fmwkout</code> makes the coding look nicer.
#define slibout sierra::Slib::theDiagWriter()

/// Macro <code>SLIB_TRACE_ENABLED</code> enables the traceback and tracing when defined.
#define SLIB_TRACE_ENABLED

#ifdef SLIB_TRACE_ENABLED
typedef Diag::Tracespec Tracespec;
typedef Diag::Traceback Traceback;

class Trace : public Diag::Trace
{
public:
  explicit Trace(const char *message)
    : Diag::Trace(slibout, message)
  {}
};
#else
typedef Diag::Tracespec Tracespec;
typedef Diag::Tracespec Traceback;
typedef Diag::Tracespec Trace;
#endif

} // namespace Slib

namespace Diag {
using stk::diag::push;
using stk::diag::pop;
using stk::diag::dendl;
} // namespace Diag

} // namespace sierra

#endif // STK_UTIL_DIAG_SlibDiagWriter_h
