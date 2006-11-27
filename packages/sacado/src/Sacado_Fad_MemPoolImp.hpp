// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include <new>

inline
Sacado::Fad::MemPool::MemPool(unsigned int elem_size, unsigned int n_elem,
			      unsigned int pre_alloc) :
  esize(elem_size < sizeof(Link) ? sizeof(Link) : elem_size),
  n(n_elem),
  csize(esize*n+sizeof(Chunk)),
  chunks(NULL),
  head(NULL),
  num_chunks(0)
{
  // Pre allocate chunks if required
  if (pre_alloc) 
    for (unsigned int i=0; i<pre_alloc; i++)
      grow();
}

inline 
Sacado::Fad::MemPool::~MemPool()
{
  Chunk *n = chunks;
  while (n != NULL) {
    Chunk *p = n;
    n = n->next;
    delete p;
  }
}

inline void*
Sacado::Fad::MemPool::alloc()
{
  if (head == NULL)
    grow();
  Link *p = head;
  head = p->next;

  return p;
}

inline void 
Sacado::Fad::MemPool::free(void *b)
{
  if (b == NULL)
    return;
  Link *p = static_cast<Link*>(b);
  p->next = head;
  head = p;
}

inline void
Sacado::Fad::MemPool::grow()
{
  // Create a new chunk
  void *p = operator new(csize);
  Chunk *c = static_cast<Chunk*>(p);
  c->mem = static_cast<char*>(p)+sizeof(Chunk);
  c->next = chunks;
  chunks = c;
  ++num_chunks;

  // Initialize each element in a chunk
  char *start = c->mem;
  char *last = &start[(n-1)*esize];
  
  for (char *q = start; q<last; q += esize)
    reinterpret_cast<Link*>(q)->next = reinterpret_cast<Link*>(q+esize);
  reinterpret_cast<Link*>(last)->next = NULL;
  head = reinterpret_cast<Link*>(start);
}
