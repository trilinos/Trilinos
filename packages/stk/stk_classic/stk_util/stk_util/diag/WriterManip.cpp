/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/diag/WriterManip.hpp>

namespace stk_classic {
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
} // namespace stk_classic
