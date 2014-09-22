/*--------------------------------------------------------------------*/
/*    Copyright (c) 2013, Sandia Corporation.
/*    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*    the U.S. Governement retains certain rights in this software.
/*    
/*    Redistribution and use in source and binary forms, with or without
/*    modification, are permitted provided that the following conditions are
/*    met:
/*    
/*        * Redistributions of source code must retain the above copyright
/*          notice, this list of conditions and the following disclaimer.
/*    
/*        * Redistributions in binary form must reproduce the above
/*          copyright notice, this list of conditions and the following
/*          disclaimer in the documentation and/or other materials provided
/*          with the distribution.
/*    
/*        * Neither the name of Sandia Corporation nor the names of its
/*          contributors may be used to endorse or promote products derived
/*          from this software without specific prior written permission.
/*    
/*    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*    
/*--------------------------------------------------------------------*/

#include <stk_util/util/Pool.hpp>

namespace stk {
namespace util {

Pool::Pool(unsigned int sz)
 : chunks(NULL),
   esize(sz<sizeof(Link) ? sizeof(Link) : sz),
   head(NULL)
{
}

Pool::~Pool()
{
  //free all chunks
  Chunk* n = chunks;
  while(n) {
    Chunk* p = n;
    n = n->next;
    delete p;
  }
}

void
Pool::grow()
{
  //allocate new chunk, organize it as a linked list of elements of size 'esize'
  Chunk* n = new Chunk;
  n->next = chunks;
  chunks = n;

  const int nelem = Chunk::size/esize;
  char* start = n->mem;
  char* last = &start[ (nelem-1)*esize ];
  for(char* p=start; p<last; p+=esize) {
    reinterpret_cast<Link*>(p)->next = reinterpret_cast<Link*>(p+esize);
  }
  reinterpret_cast<Link*>(last)->next = NULL;
  head = reinterpret_cast<Link*>(start);
}

}//namespace util
}//namespace stk

