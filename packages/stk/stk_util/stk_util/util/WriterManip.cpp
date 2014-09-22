/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/util/WriterManip.hpp>
#include <iosfwd>                       // for ostream
#include "stk_util/util/Writer.hpp"     // for Writer

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

Writer &
operator<<(Writer &dout, _resetiosflags reset_flags) {
  if (dout.shouldPrint())
    dout.getStream().setf(std::ios_base::fmtflags(0), reset_flags.m_flags);
  return dout;
}

Writer &
operator<<(Writer &dout, _setiosflags set_flags) {
  if (dout.shouldPrint())
    dout.getStream().setf(set_flags.m_flags);
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
