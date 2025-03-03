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


#ifndef STK_STK_UTIL_STK_UTIL_UTIL_DEBUGSTREAM_HPP_
#define STK_STK_UTIL_STK_UTIL_UTIL_DEBUGSTREAM_HPP_

#include <ostream>

namespace stk
{
namespace mesh
{
class NullStream : public std::ostream {
    class NullBuffer : public std::streambuf {
    public:
        int overflow( int c ) override { return c; }
    } m_nb;
public:
    NullStream() : std::ostream( &m_nb ) {}
};

class DebugStream
{
public:
  DebugStream()
  : myRank(-1)
  , outputStream(nullptr)
  { }

  DebugStream(int rank)
  : myRank(rank)
  , outputStream(nullptr)
  { }

  DebugStream(int rank, std::ostream* stream)
  : myRank(rank)
  , outputStream(stream)
  { }

  void set_output_stream(std::ostream& ostrm) { outputStream = &ostrm; }

  std::ostream& debug(bool prefixMsg = true) {
    if(outputStream != nullptr) {
      if(prefixMsg) {
        *outputStream << "P[" << myRank << "] ";
      }
      return *outputStream;
    } else {
      return ns;
    }
  }

  std::ostream& debug_p0(bool prefixMsg = true) {
    if(myRank == 0) {
      return debug(prefixMsg);
    }
    return ns;
  }

private:
  NullStream ns;
  int myRank;
  std::ostream* outputStream = nullptr;
};


}
}

#endif
