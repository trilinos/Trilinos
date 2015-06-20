// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#ifndef _stk_util_util_Pool_hpp_
#define _stk_util_util_Pool_hpp_

#include <cstdlib>

#ifndef STK_UTIL_POOL_ALLOC_CHUNK_SIZE_K
#define STK_UTIL_POOL_ALLOC_CHUNK_SIZE_K 512
#endif

//The macro STK_UTIL_POOL_ALLOC_CHUNK_SIZE_K determines the number
//of kilobytes that each internally-allocated chunk of memory will
//occupy. The Pool object will then dispense "sub-chunks"
//of memory having size determined by the argument to the class
//constructor.

namespace stk {
namespace util {
class Pool {
 public:
  Pool(unsigned int nbytes); // nbytes is the size of elements
  ~Pool();

  void* alloc(); //allocate one element
  void free(void* b); //put an element back into the pool
  struct Link { Link* next; };

 private:
  struct Chunk {
    //Stroustrup's comment:
    //slightly less than specified K so that a chunk will fit in
    //allocation area first to get stringent alignment
    enum { size = STK_UTIL_POOL_ALLOC_CHUNK_SIZE_K*1024-16 };
    char mem[size];
    Chunk* next;
  };

  Chunk* chunks;
  const unsigned int esize;
  Link* head;

  Pool(const Pool&);//private copy constructor
  Pool& operator=(const Pool&);//private assignment operator
  void grow(); //make pool larger
};

inline void* Pool::alloc()
{
  if (head == NULL) {
    grow();
  }
  Link* p = head; //return first element
  head = p->next;
  return p;
}

inline void Pool::free(void* b)
{
  Link* p = static_cast<Link*>(b);
  p->next = head; //put b back as first element
  head = p;
}

}//namespace util
}//namespace stk

#endif

