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

#include "stk_util/util/WriterManip.hpp"
#include "stk_util/util/Writer.hpp"  // for Writer
#include <iosfwd>                    // for ostream

namespace stk {
namespace diag {

Writer &operator<<(Writer &dout, _setw set_width) {
  if (dout.shouldPrint())
    dout.getStream().width(set_width.m_width);
  return dout;
}

Writer &
operator<<(Writer &dout, _setprecision set_precision) {
  if (dout.shouldPrint())
    dout.getStream().precision(set_precision.m_precision);
  return dout;
}

Writer &
operator<<(Writer &dout, _setfill set_fill) {
  if (dout.shouldPrint())
    dout.getStream().fill(set_fill.m_fill);
  return dout;
}

bool is_floatfield_hexfloat_or_defaultfloat(std::ios_base::fmtflags flags)
{
    bool neither_floatfield_bit_set = (std::ios_base::fmtflags(0) == (std::ios_base::floatfield & flags));
    bool both_floatfield_bits_set = ((std::ios_base::fixed & flags) && (std::ios_base::scientific & flags));

    return neither_floatfield_bit_set || both_floatfield_bits_set;
}

void reset_floatfield_flags_assuming_hexfloat_state(Writer &dout, std::ios_base::fmtflags reset_flags)
{
    if ((reset_flags & std::ios_base::scientific) && !(reset_flags & std::ios_base::fixed))
        fixed(dout);
    else if ((reset_flags & std::ios_base::fixed) && !(reset_flags & std::ios_base::scientific))
        scientific(dout);
}

Writer &
operator<<(Writer &dout, _resetiosflags reset_flags)
{
  if (dout.shouldPrint())
  {
    bool stream_floatfield_was_hexfloat_or_defaultfloat =
            is_floatfield_hexfloat_or_defaultfloat(dout.getStream().flags());

    dout.getStream().unsetf(reset_flags.m_flags);

    if (stream_floatfield_was_hexfloat_or_defaultfloat)
        reset_floatfield_flags_assuming_hexfloat_state(dout, reset_flags.m_flags);
  }
  return dout;
}

Writer &
operator<<(Writer &dout, _setiosflags set_flags) {
  if (dout.shouldPrint())
  {
    dout.getStream().setf(set_flags.m_flags);

    // As long as the compilers (and libstdc++ versions) we support have inconsistent
    // support for the ios_base floatfield bits, automatically do hexfloat-to-defaultfloat
    // transition.
    std::ios_base::fmtflags newFlags = dout.getStream().flags();
    if ((newFlags & std::ios_base::fixed) && (newFlags & std::ios_base::scientific))
        dout.getStream().unsetf(std::ios_base::floatfield);
  }
  return dout;
}

Writer &
fixed(
  Writer &      dout)
{
  if (dout.shouldPrint())
    dout.getStream().setf(std::ios_base::fixed, std::ios_base::floatfield);
  return dout;
}

Writer &
scientific(
  Writer &      dout)
{
  if (dout.shouldPrint())
    dout.getStream().setf(std::ios_base::scientific, std::ios_base::floatfield);
  return dout;
}

Writer &
dec(
  Writer &      dout)
{
  if (dout.shouldPrint())
    dout.getStream().setf(std::ios_base::dec, std::ios_base::basefield);
  return dout;
}

Writer &
hex(
  Writer &      dout)
{
  if (dout.shouldPrint())
    dout.getStream().setf(std::ios_base::hex, std::ios_base::basefield);
  return dout;
}

Writer &
oct(
  Writer &      dout)
{
  if (dout.shouldPrint())
    dout.getStream().setf(std::ios_base::oct, std::ios_base::basefield);
  return dout;
}

} // namespace diag
} // namespace stk
