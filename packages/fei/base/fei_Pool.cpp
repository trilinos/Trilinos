/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include "fei_macros.hpp"
#include "fei_Pool.hpp"

fei_Pool::fei_Pool(unsigned int sz)
 : chunks(NULL),
   esize(sz<sizeof(Link) ? sizeof(Link) : sz),
   head(NULL)
{
}

fei_Pool::~fei_Pool()
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
fei_Pool::grow()
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

