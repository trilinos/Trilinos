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


#ifndef STK_UTIL_DIAG_WRITER_FWD_HPP
#define STK_UTIL_DIAG_WRITER_FWD_HPP

namespace stk {

///
/// @addtogroup diag_writer_detail
/// @{
///

/**
 * @brief Enumeration <code>LogMask</code> enumerates global bit assignments.
 *
 */
enum LogMask {
  LOG_ALWAYS    = 0x00000000,
  LOG_TRACE    = 0x00000001,
  LOG_TRACE_STATS  = 0x00000002,
  LOG_TRACE_SUB_CALLS  = 0x00000004,
  LOG_MEMBERS    = 0x00000008,

  LOG_STREAM_COMMON  = LOG_TRACE | LOG_TRACE_STATS,
  LOG_PARAMETERS	= 0x00000100
};

///
/// @}
///

namespace diag {

class Writer;

typedef unsigned long PrintMask;

} // namespace diag
} // namespace stk

namespace sierra {

enum LogMask {
  LOG_STREAM_COMMON  = stk::LOG_STREAM_COMMON,
  LOG_PARAMETERS  = stk::LOG_PARAMETERS
};

namespace Diag {

typedef stk::diag::Writer Writer;

class PrintTable;

typedef stk::diag::PrintMask PrintMask;

} // namespace Diag

namespace Slib {
enum {
  LOG_ALWAYS		= stk::LOG_ALWAYS,
  LOG_TRACE		= stk::LOG_TRACE,
  LOG_TRACE_STATS	= stk::LOG_TRACE_STATS,
  LOG_TRACE_SUB_CALLS	= stk::LOG_TRACE_SUB_CALLS,
  LOG_MEMBERS		= stk::LOG_MEMBERS,

  LOG_RESOURCE		= 0x00000100,
  LOG_PLUGIN		= 0x00000200,
  LOG_GLOBAL_VARIABLE	= 0x00000400,
  LOG_MEMORY            = 0x00000800
};

} // namespace Slib

} // namespace sierra


#endif // STK_UTIL_DIAG_WRITER_FWD_HPP
