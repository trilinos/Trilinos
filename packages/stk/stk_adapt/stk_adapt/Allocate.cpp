/*
 * Allocate.cpp
 *
 *  Created on: Jan 22, 2011
 *      Author: kdcopps
 */

#ifdef STK_ADAPT_REPLACE_NEW_DELETE

#ifdef STK_HAVE_TBB
#include <tbb/scalable_allocator.h>

#pragma GCC visibility push(hidden)

extern "C++"
{
    void* operator new(std::size_t) throw (std::bad_alloc);
    void* operator new[](std::size_t) throw (std::bad_alloc);
    void* operator new(std::size_t, const std::nothrow_t&) throw();
    void* operator new[](std::size_t, const std::nothrow_t&) throw();
    void operator delete(void*) throw();
    void operator delete[](void*) throw();
    void operator delete(void*, const std::nothrow_t&) throw();
    void operator delete[](void*, const std::nothrow_t&) throw();
} // extern "C++"

#pragma GCC visibility pop


// See Intel Threading Building Blocks, James Reinders (2007).
//     chapter 11.11.1. Replacing new and delete
//
// No retry loop because we assume that scalable_malloc does
// all it takes to allocate the memory, so calling it repeatedly
// will not improve the situation at all
//
// No use of std::new_handler because it cannot be done in portable
// and thread-safe way (see sidebar)
//
// We throw std::bad_alloc() when scalable_malloc returns NULL
//(we return NULL if it is a no-throw implementation)


void* operator new (size_t size) throw (std::bad_alloc) {
    if (size == 0) size = 1;
    if (void* ptr = scalable_malloc (size))
        return ptr;
    throw std::bad_alloc ();
}

void* operator new[] (size_t size) throw (std::bad_alloc) {
    return operator new (size);
}

void* operator new (size_t size, const std::nothrow_t&) throw () {
    if (size == 0) size = 1;
    if (void* ptr = scalable_malloc (size))
        return ptr;
    return NULL;
}

void* operator new[] (size_t size, const std::nothrow_t&) throw () {
    return operator new (size, std::nothrow);
}

void operator delete (void* ptr) throw () {
    if (ptr != 0) scalable_free (ptr);
}

void operator delete[] (void* ptr) throw () {
    operator delete (ptr);
}

void operator delete (void* ptr, const std::nothrow_t&) throw () {
    if (ptr != 0) scalable_free (ptr);
}

void operator delete[] (void* ptr, const std::nothrow_t&) throw () {
    operator delete (ptr, std::nothrow);
}

#endif // STK_HAVE_TBB

#endif // STK_ADAPT_REPLACE_NEW_DELETE
