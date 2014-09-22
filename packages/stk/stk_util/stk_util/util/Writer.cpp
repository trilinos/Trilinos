// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_util/util/Writer.hpp>
#include <stk_util/util/IndentStreambuf.hpp>  // for pop, push
#include <string>                       // for string
#include "stk_util/util/Writer_fwd.hpp"  // for PrintMask


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

