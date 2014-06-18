#ifndef SAMBA_SAMBA_UTILITY_DEBUG_MESSAGE_HPP
#define SAMBA_SAMBA_UTILITY_DEBUG_MESSAGE_HPP

#include <ostream>

#ifndef NDEBUG //DEBUG -- create message

#include <sstream>
#include <stdexcept>

#include <boost/assert.hpp>

namespace boost {

// Provide our own assertion handlers to BOOST_ASSERT/BOOST_ASSERT_MSG so that we can
// get exceptions instead of aborts (helps testability).

inline
void assertion_failed(char const * expr, char const * function, char const * file, long line)
{
  std::ostringstream out;
  out << "Expression: '" << expr << "' FAILED\n"
      << "In function '" << function << "' in file '" << file << "' at line " << line << "\n";
  throw std::runtime_error(out.str());
}

inline
void assertion_failed_msg(char const * expr, char const * msg, char const * function, char const * file, long line)
{
  std::ostringstream out;
  out << "Expression: '" << expr << "' FAILED\n"
      << "In function '" << function << "' in file '" << file << "' at line " << line << "\n"
      << "Message: " << msg << "\n";
  throw std::runtime_error(out.str());
}

} // boost

namespace samba {

// TODO - Change name to dout
// TODO - Unit test to check that doesn't clobber on rethrow
// TODO - Unit test to check that it does get cleared normally

class debug_message
{
public:
  debug_message()
    : m_out()
  {}

  template <typename T>
  debug_message(T const& t)
    : m_out(t)
  {}

  const char * operator()() const
  {
    static std::string g;
    g = m_out.str();
    return g.c_str();
  }

  operator const char *() const
  {
    static std::string g;
    g = m_out.str();
    return g.c_str();
  }

  template <typename T>
  debug_message & operator<<(T const& t)
  { m_out << t; return *this; }

  debug_message & operator<<(debug_message const& msg)
  { m_out << msg.m_out.str(); return *this; }

private:
  std::ostringstream m_out;
};

inline std::ostream& operator << (std::ostream & out, debug_message const& msg)
{ return out << msg(); }

} //namespace samba

#else //RELEASE -- don't create message

namespace samba {

class debug_message
{
  public:
    debug_message()
    {}

    template <typename T>
    debug_message(T const&)
    {}

    const char * operator()() const
    { static const char * str = ""; return str; }

    operator const char *() const
    { static const char * str = ""; return str; }

    template <typename T>
    debug_message & operator<<(T const&)
    { return *this; }
};

inline std::ostream& operator << (std::ostream & out, debug_message const&)
{ return out; }

} //namespace samba

namespace boost {

inline
void assertion_failed(char const * expr, char const * function, char const * file, long line)
{}

inline
void assertion_failed_msg(char const * expr, char const * msg, char const * function, char const * file, long line)
{}

} // boost

#endif //NDEBUG

#endif //SAMBA_SAMBA_UTILITY_DEBUG_MESSAGE_HPP
