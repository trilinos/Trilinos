/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/diag/WriterExt.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/environment/Demangle.hpp>

namespace stk {
namespace diag {

Writer &
operator<<(
  Writer &                      dout,
  const std::type_info &        t)
{
  if (dout.shouldPrint())
    dout << stk::demangle(t.name());
  return dout;
}

} // namespace diag
} // namespace stk
