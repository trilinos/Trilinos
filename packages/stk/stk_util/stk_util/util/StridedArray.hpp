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

#ifndef STK_UTIL_UTIL_STRIDEDARRAY_HPP
#define STK_UTIL_UTIL_STRIDEDARRAY_HPP

#include <stk_util/stk_config.h>

namespace stk
{
namespace util
{
constexpr unsigned defaultStride = 512;

template <typename T>
class StridedArray
{
public:
  STK_FUNCTION 
  StridedArray() : dataPointer(nullptr), length(0), stride(defaultStride)
  {
  }
  STK_FUNCTION
  StridedArray(T* e, unsigned n, int stride_in=defaultStride) : dataPointer(e), length(n), stride(stride_in)
  {
  }
  STK_FUNCTION
  T operator[](unsigned i) const
  {
#ifdef __CUDA_ARCH__
    return dataPointer[stride*i];
#else
    return dataPointer[i];
#endif
  }
  STK_FUNCTION
  unsigned size() const
  { 
    return length;
  }
private:
  T* dataPointer;
  unsigned length;
  int stride;
};

} //namespace util
} //namespace stk

#endif
