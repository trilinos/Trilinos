/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
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
  typedef std::set<std::ostream *> StreamSet;
  typedef std::map<std::ostream *, int> StreamErrorMap;
  
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
   *
   * @param sb			a <b>std::streambuf</b> pointer to the output strema buffer to add.
   *
   */
  void add(std::ostream *os) {
    m_destinations.insert(os);
  }

  /**
   * @brief Member function <b>remove</b> removes the specified destination output stream buffer. 
   *
   * @param sb			a <b>std::streambuf</b> pointer to the output strema buffer to
   *                            remove. 
   *
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
