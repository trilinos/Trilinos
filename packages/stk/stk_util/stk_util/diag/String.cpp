/*
// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
 */

#include "stk_util/diag/String.hpp"
#include <cctype>    // for isspace
#include <iostream>  // for operator<<

//----------------------------------------------------------------------

namespace sierra {

namespace {

static constexpr char arraylower_t[] = {
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B,
    0x0C, 0x0D, 0x0E, 0x0F, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
    0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F, 0x20, 0x21, 0x22, 0x23,
    0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, 0x2D, 0x2E, 0x2F,
    0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B,
    0x3C, 0x3D, 0x3E, 0x3F, 0x40, 'a',  'b',  'c',  'd',  'e',  'f',  'g',
    'h',  'i',  'j',  'k',  'l',  'm',  'n',  'o',  'p',  'q',  'r',  's',
    't',  'u',  'v',  'w',  'x',  'y',  'z',  0x5B, 0x5C, 0x5D, 0x5E, 0x5F,
    0x60, 'a',  'b',  'c',  'd',  'e',  'f',  'g',  'h',  'i',  'j',  'k',
    'l',  'm',  'n',  'o',  'p',  'q',  'r',  's',  't',  'u',  'v',  'w',
    'x',  'y',  'z',  0x7B, 0x7C, 0x7D, 0x7E, 0x7F};

constexpr int arraylower(unsigned char c) noexcept {
  c = c > 127 ? 127 : c;
  return arraylower_t[c];
}


int to_label( int c )
{
  return isspace(c) ? '_' : c;
}
} // namespace <unnamed>


int char_simple_traits::compare( const char * c1 , const char * c2, size_t len ) noexcept
{
  if (len == 0) return 0;
  size_t i = 0;
  while (i < len - 1 && arraylower(c1[i]) == arraylower(c2[i]))
  {
    ++i;
  }
  return arraylower(c1[i]) - arraylower(c2[i]);
}

int char_label_traits::compare( const char * c1 , const char * c2, size_t len) noexcept
{
  if (len == 0) return 0;
  size_t i = 0;
  while (i < len - 1 && arraylower(to_label(c1[i])) == arraylower(to_label(c2[i])))
  {
    ++i;
  }
  return arraylower(to_label(c1[i])) - arraylower(to_label(c2[i]));
}

void char_label_traits::convert( char * c , size_t n )
{
  if(c==nullptr) return;
  for ( char * const e = c + n ; c != e ; ++c )
    *c = to_label(*c);
}

} // namespace sierra

//----------------------------------------------------------------------

namespace sierra {

std::ostream &
operator<<( std::ostream & os, const sierra::String & s)
{ return os << s.c_str(); }

std::istream &
operator>>( std::istream & is, sierra::String &	s )
{ std::string tmp; is >> tmp; s.assign(tmp); return is; }

std::ostream &
operator<<( std::ostream & os, const sierra::Identifier &s)
{ return os << s.c_str(); }

std::istream &
operator>>( std::istream & is, sierra::Identifier &s )
{ std::string tmp; is >> tmp; s.assign(tmp); return is; }

} // namespace sierra
