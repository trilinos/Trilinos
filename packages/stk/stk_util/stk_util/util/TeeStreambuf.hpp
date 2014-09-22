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

#ifndef STK_UTIL_UTIL_TEESTREAMBUF_HPP
#define STK_UTIL_UTIL_TEESTREAMBUF_HPP

#include <string>
#include <streambuf>
#include <ostream>
#include <set>
#include <map>

namespace stk {

namespace {
struct OStreamPointerLess
{
    inline bool operator()(const std::ostream *lhs, const std::ostream *rhs) const
    {
        return lhs < rhs;
    }
};
}
/**
 * @brief Class <b>basic_tee_streambuf</b> maintains a list of destination output stream buffers to
 * send written characters to.  Many destination output stream buffers may be added.  For each
 * character written to this stream buffer, the same character is written to each destination stream
 * buffer. 
 *
 */
template<class Ch, class Tr = std::char_traits<Ch> >
class basic_tee_streambuf : public std::basic_streambuf<Ch, Tr>
{
  typedef std::set<std::ostream *, OStreamPointerLess> StreamSet;
  typedef std::map<std::ostream *, int, OStreamPointerLess> StreamErrorMap;
  
public:
  /**
   * Creates a new <b>basic_tee_streambuf</b> instance.
   *
   */
  basic_tee_streambuf()
  {}

  /**
   * Creates a new <b>basic_tee_streambuf</b> instance and adds the specified destination output
   * stream buffer. 
   *
   */
  explicit basic_tee_streambuf(std::basic_ostream<Ch, Tr> *os) {
    add(os);
  }

  /**
   * Destroys a <b>basic_tee_streambuf</b> instance.
   *
   */
  virtual ~basic_tee_streambuf()
  {}

  /**
   * @brief Member function <b>eof</b> returns the current end-of-file status.
   *
   * @return			an <b>int</b> value of the current end-of-file status.
   */
  int eof() {
    return std::basic_streambuf<Ch, Tr>::traits_type::eof();
  }

  /**
   * @brief Member function <b>add</b> adds the specified destination output stream buffer.
   */
  void add(std::ostream *os) {
    m_destinations.insert(os);
  }

  /**
   * @brief Member function <b>remove</b> removes the specified destination output stream buffer. 
   */
  void remove(std::ostream *os) {
    m_destinations.erase(os);
  }

  /**
   * @brief Member function <b>clear</b> removes are destination output stream buffers.
   *
   */
  void clear() {
    m_destinations.clear();
  }

private:
  /**
   * @brief Member function <b>sync</b>  syncs the destination output stream buffer. 
   *
   * @return			an <b>int</b> value of 1 if successful.
   */
  virtual int sync() {
    if (m_destinations.empty())
      return 1;

    StreamErrorMap return_code;

    for (StreamSet::const_iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) {
      if ((*it)->rdbuf() != this) {
        int ret = (*it)->rdbuf()->pubsync();
        return_code[*it] = ret;
      }
    }

    // Remove streambufs with errors
    for (StreamSet::iterator it = m_destinations.begin(); it != m_destinations.end(); ++it)
      if (return_code[*it] == eof())
        m_destinations.erase(it);

    if (m_destinations.empty())
      return 1;

    return 1;
  }

  /**
   * @brief Member function <b>overflow</b> writes the specified character to all the destination
   * outptu strema buffers.
   *
   * @param c			an <b>int</b> const value of the character to write.
   *
   * @return			an <b>int</b> value of the character written.
   */
  virtual typename std::basic_streambuf<Ch, Tr>::int_type overflow(const int c) {
    if (m_destinations.empty())
      return 1;

    StreamErrorMap return_code;

    for (StreamSet::const_iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) {
      int ret = (*it)->rdbuf()->sputc(c);
      return_code[*it] = ret;
    }

    // Remove streambufs with errors
    for (StreamSet::iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) 
      if (return_code[*it] == eof())
        m_destinations.erase(it);

    if (m_destinations.empty())
      return 1;

    return 1;
  }

  /**
   * @brief Member function <b>xsputn</b> writes the specified characters to all the destination
   *
   * @param buffer		a <b>char</b> const pointer to the character string to write.
   *
   * @param n			a <b>std::streamsize</b> value of the number of characters to write.
   *
   * @return			a <b>std::streamsize</b> value of the number of characters written.
   */
  virtual std::streamsize xsputn(char const *buffer, std::streamsize n) {
    if (m_destinations.empty())
      return n;

    StreamErrorMap return_code;
    
    for (StreamSet::const_iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) {
      std::ostream *os = (*it);
      int ret = os->rdbuf()->sputn(buffer,n);
      return_code[*it] = ret;
    }

    // Remove ostreams with errors
    for (StreamSet::iterator it = m_destinations.begin(); it != m_destinations.end(); ++it) {
      if (return_code[*it] < 0) {
        m_destinations.erase(it);
      }
    }

    if (m_destinations.empty())
      return n;

    return n;
  }

private:
  StreamSet             m_destinations;    ///< Destination output stream buffers to write to
};

typedef stk::basic_tee_streambuf<char> tee_streambuf;

} // namespace stk

#endif // STK_UTIL_UTIL_TEESTREAMBUF_HPP
