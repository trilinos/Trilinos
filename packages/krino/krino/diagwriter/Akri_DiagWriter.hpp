// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_DiagWriter_h
#define Akri_DiagWriter_h

#include <stk_util/environment/Trace.hpp>
#include <stk_util/util/Writer.hpp>
#include <stk_util/diag/WriterOStream.hpp>
#include <stk_util/diag/WriterParser.hpp>

#include <Akri_DiagWriter_fwd.hpp>

namespace krino
{
  class DiagWriterParser : public stk::diag::WriterParser
  {
  public:
    DiagWriterParser();
  };

  stk::diag::Writer &theDiagWriter();
  DiagWriterParser &theDiagWriterParser();

#define krinolog krino::theDiagWriter()

#ifdef KRINO_TRACE_ENABLED
  typedef stk::diag::Tracespec Tracespec;
  typedef stk::diag::Traceback Traceback;

  class Trace : public stk::diag::Trace
  {
  public:
    Trace(const char *message)
      : stk::diag::Trace(krinolog, message)
    {}
  };

#elif defined(KRINO_TRACE_TRACEBACK)
  typedef stk::diag::Tracespec Tracespec;
  typedef stk::diag::Traceback Traceback;
  typedef stk::diag::Traceback Trace;

#else
  typedef stk::diag::Tracespec Tracespec;
  typedef stk::diag::Tracespec Traceback;
  typedef stk::diag::Tracespec Trace;
#endif

#define ThrowRuntimeError(message)                                \
{                                                                 \
  std::ostringstream internal_throw_runtime_message;              \
  internal_throw_runtime_message << message;                      \
  throw std::runtime_error(internal_throw_runtime_message.str()); \
}
} // namespace krino

#endif // Akri_DiagWriter_h
