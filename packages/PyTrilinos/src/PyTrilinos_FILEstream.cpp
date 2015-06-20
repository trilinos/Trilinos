// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "PyTrilinos_FILEstream.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <cstdio>

using std::size_t;

namespace PyTrilinos
{

FILEstream::FILEstream(FILE   *fptr,
		       size_t buff_sz,
		       size_t put_back) :
  fptr_(fptr),
  put_back_(std::max(put_back, size_t(1))),
  buffer_(std::max(buff_sz, put_back_) + put_back_)
{
  char *beg = &buffer_.front();
  char *end = beg + buffer_.size();
  setg(end, end, end);    // Set the get buffer
  setp(beg, end-1);       // Set the put buffer
}

std::streambuf::int_type FILEstream::underflow()
{
  if (gptr() < egptr())        // buffer not exhausted
    return traits_type::to_int_type(*gptr());

  char *base  = &buffer_.front();
  char *start = base;

  if (eback() == base)         // true if not the first fill
  {
    // Make arrangements for putback characters
    std::memmove(base, egptr() - put_back_, put_back_);
    start += put_back_;
  }

  size_t n = std::fread(start, 1, buffer_.size() - (start - base), fptr_);
  if (n == 0)
    return traits_type::eof();

  setg(base, start, start+n);
  return traits_type::to_int_type(*gptr());
}

std::streambuf::int_type FILEstream::overflow(std::streambuf::int_type ch)
{
  sync();
  if (ch != traits_type::eof())
  {
    *pptr() = ch;
    pbump(1);
  }
  return ch;
}

std::streambuf::int_type FILEstream::sync()
{
  char   *beg  = &buffer_.front();
  char   *end  = beg + buffer_.size();
  size_t count = pptr() - beg;
  size_t n     = std::fwrite(beg, 1, count, fptr_);
  setp(beg, end-1);
  if (n != count)
    return -1;
  return 0;
}

}  // Namespace PyTrilinos
