// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_STREAM_HPP
#define ROL_STREAM_HPP

#include <ostream>
#include <string>
#include "ROL_Ptr.hpp"

/** \file  ROL_Stream.hpp
    \brief Defines a no-output stream class ROL::NullStream and a function
           makeStreamPtr which either wraps a reference to a stream object
           or returns a pointer to a NullStream depending on the value of
           the argument noSuppressOutput
           
*/

namespace ROL {

namespace details {

using namespace std;

template<typename _CharT, typename _Traits>
class basic_nullstream : virtual public basic_ostream<_CharT, _Traits> {
public:
  explicit basic_nullstream() : basic_ostream<_CharT, _Traits>(NULL) {}
}; 

using nullstream = basic_nullstream<char, char_traits<char>>;

inline
Ptr<ostream> makeStreamPtr( ostream& os, bool noSuppressOutput=true ) {
  Ptr<ostream> retstream;
  if( noSuppressOutput ) retstream = makePtrFromRef<ostream>(os);
  else retstream = makePtr<nullstream>();
  return retstream; // noSuppressOutput ? makePtrFromRef( os ) : makePtr<nullstream>();
}

inline
Ptr<ostream> makeStreamPtr( Ptr<ostream> os, bool noSuppressOutput=true ) {
  Ptr<ostream> retstream;
  if( noSuppressOutput ) retstream = os;
  else retstream = makePtr<nullstream>();
  return retstream; // noSuppressOutput ? makePtrFromRef( os ) : makePtr<nullstream>();
//  return noSuppressOutput ? os : makePtr<nullstream>();
}

} // details

using details::nullstream;
using details::makeStreamPtr;

} // namespace ROL


#endif 
