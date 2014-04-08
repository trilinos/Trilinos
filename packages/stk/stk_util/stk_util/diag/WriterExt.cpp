/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/diag/WriterExt.hpp>
#include <stk_util/diag/String.hpp>     // for Identifier, String
#include <stk_util/environment/Demangle.hpp>  // for demangle
#include <stk_util/util/Writer.hpp>     // for operator<<, Writer
#include "stk_util/parallel/MPI.hpp"    // for Loc, TempLoc



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


Writer &
operator<<(
  Writer &                      dout,
  const sierra::String &        s)
{
  if (dout.shouldPrint()) 
    dout << s.c_str();
  return dout;
}


Writer &
operator<<(
  Writer &                      dout,
  const sierra::Identifier &    s)
{
  if (dout.shouldPrint())
    dout << s.c_str();
  return dout;
}


Writer &
operator<<(
  Writer &        dout,
  const sierra::MPI::Loc<int> &      loc) 
{
  if (dout.shouldPrint()) 
    dout << loc.m_value << "@" << loc.m_loc;
  return dout;
}


Writer &
operator<<(
  Writer &        dout,
  const sierra::MPI::Loc<double> &   loc)
{
  if (dout.shouldPrint())
    dout << loc.m_value << "@" << loc.m_loc;
  return dout;
}


Writer &
operator<<(
  Writer &        dout,
  const sierra::MPI::Loc<float> &    loc)
{
  if (dout.shouldPrint())
    dout << loc.m_value << "@" << loc.m_loc;
  return dout;
}

  
Writer &
operator<<(
  Writer &        dout,
  const sierra::MPI::TempLoc &   loc)
{
  if (dout.shouldPrint())
    dout << loc.m_value << " " << loc.m_other << "@" << loc.m_loc;
  return dout;
}


} // namespace diag
} // namespace stk
