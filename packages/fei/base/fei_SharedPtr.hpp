/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_SharedPtr_hpp_
#define _fei_SharedPtr_hpp_

#include <fei_macros.hpp>
//
//fei::SharedPtr is a copy of the Sierra system's SharedPtr class, which
//was added by Kevin Copps. This class is a close copy of the boost shared_ptr.
//
//NOTE: In this copy, I've removed the member function 'swap', and the 
//std::less specialization.
//
//boost::shared_ptr now allows a second template parameter which specifies
//a deleter object. Instead of adopting that, I'm taking a lazy approach and
//simply adding a default bool argument to a constructor which specifies 
//whether the SharedPtr should delete the pointer or not.
//

// #ifdef SIERRA_NO_MEMBER_TEMPLATES
// #define FEI_NO_MEMBER_TEMPLATES
// #endif

namespace fei {

/**
 * A smart pointer with reference counted copy semantics. This class is
 * a close copy of boost::shared_ptr. See www.boost.org, and
 * www.boost.org/libs/smart_ptr/smart_ptr.htm.
 * For logistical reasons it was deemed easier to make a copy of the
 * boost::shared_ptr rather than having FEI be dependent on boost.
 * According to TR1 (C++ Technical Report 1), referenced on the boost
 * web site, shared_ptr will be migrating into the C++ standard.
 *
 * fei::SharedPtr usage notes:
 *
 * Construction:
 * - Construct fei::SharedPtr with a raw pointer obtained from new:
 *    fei::SharedPtr<MyObject> myptr(new MyObject);
 * - fei::SharedPtr also has a copy-constructor.
 *
 * Assignment:
 * - fei::SharedPtr objects can be assigned to each other. You can not
 *   assign a raw pointer to a fei::SharedPtr. Instead, use the reset
 *   member:
 *     myptr.reset(new MyObject);
 *
 * Comparison:
 * - fei::SharedPtr objects can be compared to each other, but can
 *   not be compared to raw pointers. Instead, use the get member:
 *     if (myptr.get() == NULL) {...
 *
 * The object pointed to is deleted when the last SharedPtr pointing
 * to it is destroyed or reset.
 *
 * @author Kevin Copps
 * @date 12/4/2001
 */ 
template<typename T> class SharedPtr { 
 
  public: 
 
    /** 
     * The type of the stored pointer. 
     */ 
    typedef T element_type; 
 
    /** 
     * Constructs a SharedPtr, storing a copy of p, which must have 
     * been allocated via a C++ new expression or be 0. On exit,  
     * use_count() is 1 (even if p==0; see the destructor). 
     * 
     * The only exception which may be thrown by this constructor is  
     * std::bad_alloc. If an exception is thrown, delete p is called. 
     * 
     * @param p   the pointer value recently allocated 
     */ 
    explicit SharedPtr(T* p) 
      : pointer(p)
      {
	try { // prevent leak if new throws 
	  count = new long(1); 
	} catch (...) {
	  delete p;
	  throw;
	}
      }
 
    SharedPtr(void) 
      : pointer(0) {
        count = new long(1); 
    } 


   /** 
    * Destructor. If use_count() == 1, deletes the object pointed 
    * to by the stored pointer. Otherwise, use_count() for any 
    * remaining copies is decremented by 1. Note that in C++ 
    * delete on a pointer with a value of 0 is harmless. 
    * 
    * Never throws an exception. 
    */ 
   ~SharedPtr() { dispose(); } 
 
#if !defined( FEI_NO_MEMBER_TEMPLATES ) 
  /** 
   * Constructs a SharedPtr, as if by storing a copy of the pointer 
   * stored in x. Afterwards, use_count() for all copies is 1 more  
   * than the initial x.use_count(). 
   * 
   * Never throws an exception. 
   * 
   * @param x   a shared pointer to another type 
   */ 
  template<typename Y> 
     SharedPtr(const SharedPtr<Y>& x)
    : pointer(x.pointer)
    { 
      ++*(count = x.count);
    } 
 
  /** 
   * Assignment to a shared pointer of another type. 
   * 
   * First, if use_count() == 1, deletes the object pointed to by the  
   * stored pointer. Otherwise, use_count() for any remaining copies 
   * is decremented by 1. Note that in C++  delete on a pointer with a 
   * value of 0 is harmless. 
   * 
   * Then replaces the contents of this, as if by storing a copy of  
   * the pointer stored in x. Afterwards, use_count() for all copies  
   * is 1 more than the initial x.use_count(). 
   *  
   * Never throws an exception. 
   * 
   * @param x   a shared pointer to another type 
   */ 
  template<typename Y> 
    SharedPtr& operator=(const SharedPtr<Y>& x) {  
        share(x.pointer,x.count);
        return *this; 
    } 
#endif // FEI_NO_MEMBER_TEMPLATES 
 
  /** 
   * Constructs a SharedPtr, as if by storing a copy of the pointer 
   * stored in x. Afterwards, use_count() for all copies is 1 more  
   * than the initial x.use_count(). 
   * 
   * Never throws an exception. 
   * 
   * @param x   the shared pointer to copy 
   */ 
  SharedPtr(const SharedPtr& x) 
    : pointer(x.pointer)
    {
      ++*(count = x.count);
    }
 
  /** 
   * Assignment to another shared pointer. 
   * First, if use_count() == 1, deletes the object pointed to by the  
   * stored pointer. Otherwise, use_count() for any remaining copies 
   * is decremented by 1. Note that in C++  delete on a pointer with a 
   * value of 0 is harmless. 
   * 
   * Then replaces the contents of this, as if by storing a copy of  
   * the pointer stored in x. Afterwards, use_count() for all copies  
   * is 1 more than the initial x.use_count(). 
   *  
   * Does not throw any exception. 
   *  
   * @param x  the shared pointer to copy 
   */ 
  SharedPtr& operator=(const SharedPtr& x) { 
    share(x.pointer, x.count);
    return *this; 
  } 
 
  /** 
   * Reset the pointer value of this shared pointer. 
   * First, if use_count() == 1, deletes the object pointed to by the 
   * stored pointer. Otherwise, use_count() for any remaining copies 
   * is decremented by 1. Then replaces the contents of this, as if 
   * by storing a copy of p, which must have been allocated via a C++ 
   * new expression or be 0. Afterwards, use_count() is 1 
   * (even if p==0; see ~SharedPtr). 
   *  
   * Note that in C++ delete on a pointer with a value of 0 is 
   * harmless. 
   * 
   * The only exception which may be thrown is std::bad_alloc. 
   * If an exception is thrown, delete p is called. 
   * 
   * @param p   a pointer value, or 0 if not present 
   */ 
  void reset(T* p=0) {
    if ( pointer == p ) return;
    if (--*count == 0) { 
      delete pointer;
    } 
    else { // allocate new reference counter
      try { 
	count = new long; 
      } 
      catch (...) { 
	++*count; 
	delete p; 
	throw; 
      }
    }
    *count = 1;
    pointer = p; 
  } 
 
  /** 
   * Returns a reference to the object pointed to by the stored 
   * pointer. 
   * 
   * Never throws an exception. 
   */ 
  T& operator*() const          { return *pointer; } 
 
  /** 
   * Return the stored pointer. 
   * 
   * Never throws an exception. 
   */ 
  T* operator->() const         { return pointer; } 
 
  /** 
   * Return the stored pointer. T is not required to be a complete type. 
   * 
   * Never throws an exception. 
   */ 
  T* get() const                { return pointer; } 
 
  /** 
   * Returns the number of SharedPtr's sharing ownership  
   * of the stored pointer. T is not required to be a complete type. 
   * 
   * Never throws an exception. 
   */  
  long use_count() const        { return *count; } 
 
  /** 
   * Returns use_count() == 1. 
   * T is not required to be a complete type. 
   * 
   * Never throws an exception. 
   */  
  bool unique() const           { return *count == 1; } 
 
  /** power users only */
  void share(T* xpointer, long* xcount) { 
    if (count != xcount) { 
      ++*xcount; 
      dispose(); 
      pointer = xpointer; 
      count = xcount; 
    } 
  } 
 
  /** power users only */
  void dispose() { 
    if (--*count == 0) { 
      delete pointer; 
      delete count; 
    } 
  } 
 
  // Making all members public allows member templates 
  // to work in the absence of member template friends. 
#if defined( FEI_NO_MEMBER_TEMPLATES ) || !defined( FEI_NO_MEMBER_TEMPLATES ) 
   private: 
#endif 
 
   T*     pointer; // contained pointer 
   long*  count;   // ptr to reference counter

#if !defined( FEI_NO_MEMBER_TEMPLATES ) && !defined( FEI_NO_MEMBER_TEMPLATES ) 
   template<typename Y> friend class SharedPtr; 
#endif 
 
};  // end class SharedPtr 
 
/** 
 * Equals operator for shared pointers. 
 */ 
template<typename T, typename U> 
  inline bool operator==(const SharedPtr<T>& a, const SharedPtr<U>& b) 
    { return a.get() == b.get(); } 
 
/** 
 * Not equals operator for shared pointers. 
 */ 
template<typename T, typename U> 
  inline bool operator!=(const SharedPtr<T>& a, const SharedPtr<U>& b) 
    { return a.get() != b.get(); } 

} // namespace fei

#endif // _fei_SharedPtr_hpp_

