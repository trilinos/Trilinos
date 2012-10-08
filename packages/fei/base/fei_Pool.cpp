/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

