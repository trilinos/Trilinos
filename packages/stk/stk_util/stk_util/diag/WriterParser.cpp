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

#include "stk_util/diag/WriterParser.hpp"
#include "stk_util/environment/Trace.hpp"  // for Trace
#include "stk_util/util/Writer_fwd.hpp"    // for LOG_MEMBERS, LOG_TRACE, LOG_TRACE_SUB_CALLS
#include <iomanip>                         // for operator>>, resetiosflags
#include <iostream>                        // for basic_istream, basic_istream<>::__istream_type
#include <sstream>
#include <map>                             // for _Rb_tree_const_iterator, map<>::const_iterator
#include <stdexcept>                       // for runtime_error
#include <utility>                         // for pair



namespace stk {
namespace diag {

WriterParser::WriterParser()
  : OptionMaskParser()
{
  mask("coverage", 0, "Collect and display traceable function usage coverage");
  mask("members", LOG_MEMBERS, "Display data structure members messages");
  mask("trace", LOG_TRACE, "Display execution trace");
  mask("trace-stats", LOG_TRACE_STATS, "Display execution time and memory usage during trace");
  mask("trace-down", LOG_TRACE_SUB_CALLS, "Display subsequent calls after tracing is enabled");
}


OptionMaskParser::Mask
WriterParser::parse(
  const char *          mask_string) const
{
  m_optionMask = LOG_MEMBERS;
  return OptionMaskParser::parse(mask_string);
}

} // namespace diag
} // namespace stk
