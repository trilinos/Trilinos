// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_HANDLE_HPP
#define TEUCHOS_HANDLE_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Handleable.hpp"

namespace Teuchos 
{

  /** 
   * ConstHandle is a templated handle class with strong const protection.
   *
   * In writing derived types, it is usually simplest to use the
   * TEUCHOS_CONST_HANDLE_CTORS macro to generate boilerplate constructor
   * code.
   * 
   * There are two modes of construction:
   * construction from an existing RefCountPtr,
   * \code
   * RefCountPtr<const Base> r = rcp(new Derived(blahblah));
   * ConstHandle<Base> h = r;
   * \endcode
   * and construction from a raw pointer,
   * \code
   * ConstHandle<Base> h = new Derived(blahblah);
   * \endcode
   * The second form makes the code slightly cleaner. 
   * Note that to use this second form, it is necessary that Derived implement the
   * ConstHandleable interface; this is necessary to avoid any implicit conversions
   * from raw pointers to smart pointers.
   * 
   * Note that the first form must be used whenever the object being handled has
   * been allocated on the stack.
   *
   */
  template <typename PointerType> 
  class ConstHandle : public virtual Describable
  {
  public:
    /** Construct with an existing RefCountPtr */
    explicit ConstHandle(const RefCountPtr<const PointerType>& ptr) : ptr_(ptr) {;}
    /** Construct with a raw pointer to a ConstHandleable. This will make
     * a call to rcp(), thus removing it from the user interface. */
    explicit ConstHandle(const ConstHandleable<PointerType>* ptr) : ptr_(ptr->getConstRcp()) {;}

    /** Read-only access to the underlying smart pointer. */
    const RefCountPtr<const PointerType>& ptr() const {return ptr_;}


  protected:
    /** The empty ctor will only be called by Handle ctors */
    explicit ConstHandle() : ptr_() {;}
    /** This function is needed in Handle ctors. The Handle ctors call the empty
     * ConstHandle ctor and then set the pointer in the ConstHandle with a call
     * to setRcp(). */
    void setRcp(const RefCountPtr<PointerType>& ptr) 
    {ptr_=rcp_const_cast<const PointerType>(ptr);}

    

    /** Protected non-const access to the underlying smart pointer. */
    RefCountPtr<PointerType> nonConstPtr() 
    {return rcp_const_cast<PointerType>(ptr_);}
    
  private:
    /** */
    RefCountPtr<const PointerType> ptr_;
  };


  /** 
   *  Handle is a generic templated handle class.
   *
   * In writing derived types, it is usually simplest to use the
   * TEUCHOS_HANDLE_CTORS macro to generate boilerplate constructor
   * code.
   * 
   * There are two modes of construction:
   * construction from an existing RefCountPtr,
   * \code
   * RefCountPtr<Base> r = rcp(new Derived(blahblah));
   * Handle<Base> h = r;
   * \endcode
   * and construction from a raw pointer,
   * \code
   * Handle<Base> h = new Derived(blahblah);
   * \endcode
   * The second form makes the code slightly cleaner. 
   * Note that to use this second form, it is necessary that Derived implement the
   * Handleable interface; this is necessary to avoid any implicit conversions
   * from raw pointers to smart pointers.
   * 
   * Note that the first form must be used whenever the object being handled has
   * been allocated on the stack.
   */
  template <typename PointerType> 
  class Handle : public virtual ConstHandle<PointerType>
  {
  public:
    /** Empty ctor */
    Handle()
      : ConstHandle<PointerType>() {}

    /** Construct with an existing RefCountPtr */
    explicit Handle(const RefCountPtr<PointerType>& smartPtr)
      : ConstHandle<PointerType>() 
    {
      /* We need to set the rcp in the base class */
      setRcp(smartPtr);
    }

    /** Construct with a raw pointer to a Handleable. This will make
     * a call to rcp() internally, thus removing it from the user interface. */
    explicit Handle(Handleable<PointerType>* rawPtr)
      : ConstHandle<PointerType>()
    {
      /* We need to set the rcp in the base class */
      setRcp(rawPtr->getRcp());
    }

    /** Read/write access to the underlying smart pointer. 
     *
     * \note This returns by value, not by reference. Returning by reference would
     * be <b> very </b> dangerous in this context, because it would be possible
     * to assign to the return value, e.g.,
     * \code
     * someHandle.ptr() = rcp(Blah());
     * \endcode
     * That would have the effect of desynchronizing the non-const
     * and const pointers stored by this object, with catastrophic consequences. 
     * Returning by value prevents this
     * disaster.
     */
    RefCountPtr<PointerType> ptr() {return this->nonConstPtr();}

  };

  

}

/** \def This helper macro defines boilerplate constructors 
 * for classes deriving
 * from Handle.
 *
 * If class MyHandle is a handle to a type MyType, simply 
 * put
 * \code
 * TEUCHOS_HANDLE_CTORS(MyHandle, MyType);
 * \endcode
 * in the class declaration of MyHandle and the macro will create 
 * an empty ctor, a ctor from a smart ptr, and a ctor from a raw pointer. 
 * The macro will also create appropriate doxygen for the handle ctors */
#define TEUCHOS_HANDLE_CTORS(handle, contents) \
/** Empty ctor */ \
handle() : Teuchos::Handle<contents >() {;} \
/** Construct a #handle with a raw pointer to a #contents */ \
handle(Teuchos::Handleable<contents >* rawPtr) : Teuchos::Handle<contents >(rawPtr) {;} \
/** Construct a #handle with a smart pointer to a #contents */ \
handle(const Teuchos::RefCountPtr<contents >& smartPtr) : Teuchos::Handle<contents >(smartPtr){;}

/** \def This helper macro defines boilerplate constructors 
 * for classes deriving
 * from ConstHandle.
 *
 * If class MyHandle is a const handle to a type MyType, simply 
 * put
 * \code
 * TEUCHOS_CONST_HANDLE_CTORS(MyHandle, MyType);
 * \endcode
 * in the class declaration of MyHandle and the macro will create 
 * an empty ctor, a ctor from a smart ptr, and a ctor from a raw pointer. 
 * The macro will also create appropriate doxygen for the handle ctors */
#define TEUCHOS_CONST_HANDLE_CTORS(handle, contents) \
/** Empty ctor */ \
handle() : Teuchos::ConstHandle<contents >() {;} \
/** Construct a #handle with a raw pointer to a #contents */ \
handle(const Teuchos::ConstHandleable<contents >* rawPtr) : Teuchos::ConstHandle<contents >(rawPtr) {;} \
/** Construct a #handle with a smart pointer to a #contents */ \
handle(const Teuchos::RefCountPtr<const contents >& smartPtr) : Teuchos::ConstHandle<contents >(smartPtr){;}




#endif // TEUCHOS_CONSTHANDLE_HPP
