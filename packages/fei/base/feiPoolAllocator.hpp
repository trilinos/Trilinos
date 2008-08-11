/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _feiPoolAllocator_hpp_
#define _feiPoolAllocator_hpp_

#include <fei_macros.hpp>

/** Simple Pool Allocator. Templated on the type of the objects being allocated.
When the user of this class requests the first object pointer, a pool of size
poolsize objects is allocated and a pointer to the first object in the pool is
returned. When subsequent object pointers are requested, the user is given a
pointer to the next object in the pool. A new pool or 'chunk' of memory is
allocated when necessary. Currently, none of the pools are destroyed until the
feiPoolAllocator is destroyed.
The type 'T' on which this class is templated must have a default constructor.
*/

template<class T>
class feiPoolAllocator {
 public:
  /** Constructor.
      @param poolsize Optional argument, defaults to 1024 if not provided.
      poolsize is the number of objects that will be allocated in the pool, in
      each "chunk". i.e., it is the amount by which the pool grows each time it
      is expanded.
  */
  feiPoolAllocator(int poolsize = 1024)
    : n_(poolsize),
      current_(n_),
      mempool(0)
    {
      //insist on a pool-increment-size of at least 1
      if (n_ < 1) {
	n_ = 1;
      }
    }

  /** Destructor. */
  virtual ~feiPoolAllocator()
    {
      memLink* nextLink = 0;
      if (mempool != 0) {
	delete [] mempool->memData;
	nextLink = mempool->next;
	delete mempool;
      }

      while(nextLink) {
	memLink* thisLink = nextLink;
	delete [] thisLink->memData;
	nextLink = thisLink->next;
	delete thisLink;
      }
    }

  /** Provide a pointer to one object of type T.*/
  T* alloc()
    {
      if (current_ == n_) {
	expandPool();
      }

      T* ptr = &(mempool->memData[current_]);
      current_++;
      return( ptr );
    }

  /** Provide a pointer to one or more objects of type T.
   @param num Input. Number of items to be marked as used in the pool. For
  example, if the type T is int, then by default alloc() should be assumed to
  return a pointer to a single int. But if num is specified to be greater than 1,
  then alloc() returns a pointer to 'num' ints.*/
  T* alloc(int num)
    {

      if (current_+num > n_) {
	expandPool();
      }

      T* ptr = &(mempool->memData[current_]);
      current_+=num;
      return( ptr );
    }

 private:
  feiPoolAllocator(const feiPoolAllocator<T>& src)
    : n_(0), current_(0), mempool(0) {}

  feiPoolAllocator<T>& operator=(const feiPoolAllocator<T>& src)
    { return(*this); }

  struct memLink { T* memData; memLink* next; };

  void expandPool()
    {
      memLink* newChunk = new memLink;

      newChunk->memData = new T[n_];

      newChunk->next = mempool;
      mempool = newChunk;
      current_ = 0;
    }
 
  int n_;
  int current_;

  memLink* mempool;
};

#endif // _feiPoolAllocator_hpp_
