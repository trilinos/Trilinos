#pragma once
#ifndef ROL2_STREAM_DECL_HPP
#define ROL2_STREAM_DECL_HPP

/** \file  ROL2_Stream.hpp
    \brief Defines a no-output stream class ROL::NullStream and a function
           makeStreamPtr which either wraps a reference to a stream object
           or returns a pointer to a NullStream depending on the value of
           the argument noSuppressOutput

*/

namespace ROL2 {

template<typename _CharT, typename _Traits>
class BasicNullStream : virtual public std::basic_ostream<_CharT, _Traits> {
public:
  explicit BasicNullStream() : std::basic_ostream<_CharT, _Traits>(nullptr) {}
};

using NullStream = BasicNullStream<char, std::char_traits<char>>;

inline Ptr<std::ostream> makeStreamPtr( std::ostream&     os, 
                                        bool              noSuppressOutput=true );

inline Ptr<std::ostream> makeStreamPtr( Ptr<std::ostream> os, 
                                        bool              noSuppressOutput=true );

} // namespace ROL2

#endif // ROL2_STREAM_DECL_HPP

