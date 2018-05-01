// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL_STREAM_HPP
#define ROL_STREAM_HPP

#include <ostream>
#include <streambuf>
#include <type_traits>

/** \file  ROL_Stream.hpp
    \brief Defines a no-output stream class ROL::NullStream and a function
           makeStreamPtr which either wraps a reference to a stream object
           or returns a pointer to a NullStream depending on the value of
           the argument noSuppressOutput
           
*/

namespace ROL {

namespace details {

using namespace std;

class NullBuffer : virtual public streambuf {
public:
  int overflow( int c ) override { return c; }
}; // NullBuffer

class NullStream : virtual public ostream {
private:
  NullBuffer buffer; 
public:
  NullStream( ) : ostream(&buffer) {}
}; // NullStream

inline
Ptr<ostream> makeStreamPtr( ostream& os, bool noSuppressOutput=true ) {
  Ptr<ostream> osPtr;
  if( noSuppressOutput ) osPtr = makePtrFromRef( os );
  else {
    osPtr = makePtr<NullStream>();
  }
  return osPtr;
}

} // namespace details

using details::NullStream;
using details::makeStreamPtr;

} // namespace ROL


#endif 
