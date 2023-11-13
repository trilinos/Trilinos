// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_DiagWriter_fwd_h
#define Akri_DiagWriter_fwd_h

#define KRINO_TRACE_ENABLED
// #define KRINO_TRACE_TRACEBACK

#include <stk_util/util/Writer_fwd.hpp>

namespace krino
{
  class DebugOut;

  enum
    {
      LOG_ALWAYS          = stk::LOG_ALWAYS,
      LOG_MEMBERS         = stk::LOG_MEMBERS,
      LOG_TRACE           = stk::LOG_TRACE,
      LOG_TRACE_STATS     = stk::LOG_TRACE_STATS,
      LOG_TRACE_SUB_CALLS = stk::LOG_TRACE_SUB_CALLS,

      LOG_DEBUG           = 0x00010000,
      LOG_SUBELEMENT      = 0x00020000,
      LOG_FACETS          = 0x00040000,
      LOG_PARTS           = 0x00080000
    };

} // namespace krino

#endif // Akri_DiagWriter_fwd_h
