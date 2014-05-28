/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>
#include <string>
#include <iomanip>
#include <list>

#include <stk_util/util/IndentStreambuf.hpp>
#include <stk_util/diag/Writer.hpp>

namespace stk {
namespace diag {

WriterThrowSafe::WriterThrowSafe(
  Writer &              writer)
  : m_writer(writer),
    m_depth(writer.getDepth())
{}


WriterThrowSafe::~WriterThrowSafe()
{
  m_writer.restoreDepth(m_depth);
}


Writer::Writer(
  std::streambuf *      writer_streambuf,
  PrintMask    print_mask,
  Flags      flags)
  : m_flags(flags),
    m_printMask(print_mask),
    m_lineMaskStack(),
    m_traceDepth(0),
    m_writerStream(writer_streambuf)
{}


Writer::~Writer()
{}


/**
 * @brief Member function <b>dflush</b> flushes the output stream.
 *
 * @return      a <b>Writer</b> reference to this object
 */
Writer &
Writer::dflush() {
  getStream() << std::flush;
  return *this;
}

/**
 * @brief Member function <b>dendl</b> is a manipulator which sets the output
 * stream to a new line.<P>
 *
 * The std::endl manipulator is sent to the output stream.
 *
 * @return      a <b>Writer</b> reference to this object
 */
Writer &
Writer::dendl() {
  if (shouldPrint())
    getStream() << std::endl;

  m_lineMaskStack.resetDepth();

  return *this;
}

/**
 * @brief Member function <b>push</b> is a manipulator which increases the line
 * mask depth by one.
 *
 * @return      a <b>Writer</b> reference to this object
 */
Writer &
Writer::push() {
  if (shouldPrint()) {
    m_lineMaskStack.pushDepth();
    getStream() << stk::push;
  }

  return *this;
}

/**
 * @brief Member function <b>pop</b> is a manipulator which decreases the line
 * mask depth by one, but not less than zero(0).
 *
 * @return      a <b>Writer</b> reference to this object
 */
Writer &
Writer::pop() {
  if (shouldPrint()) {
    getStream() << stk::pop;
    m_lineMaskStack.resetDepth().pop();
  }

  return *this;
}

/**
 * @brief Member function <b>pop</b> is a manipulator which decreases the line
 * mask depth by one, but not less than zero(0).
 *
 * @return      a <b>Writer</b> reference to this object
 */
Writer &
Writer::resetLineMask() {
  m_lineMaskStack.popLineMask();

  return *this;
}

/**
 * @brief Member function <b>operator<<</b> is the manipulator instantiation
 * function
 *
 * @return      a <b>Writer</b> reference to this object
 */
Writer &
Writer::operator<<(Writer& (*f)(Writer&)) {
  f(*this);
  return *this;
}

Writer &
Writer::operator<<(
  std::ios_base &       (*f)(std::ios_base&))
{
  if (shouldPrint())
    f(getStream());
  return *this;
}


Writer &
Writer::operator<<(
  std::ostream &        (*f)(std::ostream&))
{
  if (shouldPrint())
    f(getStream());

  return *this;
}


Writer &
operator<<(
  Writer &  dout,
  const void *  ptr)
{
  if (dout.shouldPrint())
    dout.getStream() << ptr;

  return dout;
}


Writer &
operator<<(
  Writer &  dout,
  const char *  c_str)
{
  if (dout.shouldPrint()) {
    std::ostream &os = dout.getStream();
    if (!c_str)
      os << "(null)";
    else
      os << c_str;
  }

  return dout;
}


Writer &
operator<<(
  Writer &              dout,
  const std::string &   s)
{
  dout << s.c_str();
  return dout;
}


Writer &
operator<<(
  Writer &  dout,
  const float & x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &          dout,
  const double &        x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &          dout,
  const long double &   x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &  dout,
  const int &   x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &          dout,
  const unsigned int &  x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &  dout,
  const long &  x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &          dout,
  const unsigned long & x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}

Writer &
operator<<(
  Writer &            dout,
  const long long &     x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}

Writer &
operator<<(
  Writer &            dout,
  const unsigned long long &  x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}

Writer &
operator<<(
  Writer &  dout,
  const short & x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


Writer &
operator<<(
  Writer &            dout,
  const unsigned short &  x)
{
  if (dout.shouldPrint())
    dout.getStream() << x;

  return dout;
}


} // namespace diag
} // namespace stk

namespace sierra {
namespace Diag {

Writer &
operator<<(
  Writer &	        dout,
  const String &        str)
{
  if (dout.shouldPrint())
    dout.getStream() << str;

  return dout;
}


Writer &
operator<<(
  Writer &                      dout,
  const sierra::Identifier &    s)
{
  if (dout.shouldPrint())
    dout.getStream() << '\'' << s << '\'';

  return dout;
}

} // namespace Diag
} // namespace sierra

