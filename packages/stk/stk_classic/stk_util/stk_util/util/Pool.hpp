/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

