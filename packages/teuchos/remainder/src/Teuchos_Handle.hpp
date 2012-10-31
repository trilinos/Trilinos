// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_HANDLE_HPP
#define TEUCHOS_HANDLE_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Handleable.hpp"

namespace Teuchos 
{

  /** \brief Templated handle class with strong const protection.
   *
   * In writing derived types, it is usually simplest to use the
   * TEUCHOS_CONST_HANDLE_CTORS macro to generate boilerplate constructor
   * code.
   * 
   * There are two modes of construction:
   * construction from an existing RCP,
   * \code
   * RCP<const Base> r = rcp(new Derived(blahblah));
   * ConstHandle<Base> h = r;
   * \endcode
   * and construction from a raw pointer,
   * \code
   * ConstHandle<Base> h = new Derived(blahblah);
   * \endcode
   * The second form makes the code slightly cleaner. Note that to use this
   * second form, it is necessary that Derived implement the ConstHandleable
   * interface; this is necessary to avoid any implicit conversions from just
   * any raw pointer to a smart pointer.
   * 
   * Note that the first form with rcp() must be used whenever the object
   * being handled has been allocated on the stack (using rcp(ptr,false) of
   * course).
   */
  template <typename PointerType> 
  class ConstHandle : public virtual Describable
  {
  public:
    /** \brief Construct with an existing RCP. */
    ConstHandle(const RCP<const PointerType>& ptr) : ptr_(ptr) {;}
    /** \brief Construct with a raw pointer to a ConstHandleable. This will make
     * a call to rcp(), thus removing that call from the user interface. */
    explicit ConstHandle(const ConstHandleable<PointerType>* ptr) : ptr_(ptr->getConstRcp()) {;}
    /** \brief Read-only access to the underlying smart pointer. */
    const RCP<const PointerType>& constPtr() const {return ptr_;}
    /** \brief Access to raw pointer */
    const PointerType * const rawPtr() {return this->constPtr().get();}
  protected:
    /** \brief The empty ctor will only be called by Handle ctors */
    explicit ConstHandle() : ptr_() {;}
    /** \brief This function is needed in Handle ctors.
     *
     * The Handle ctors call the empty ConstHandle ctor and then set the
     * pointer in the ConstHandle with a call to setRcp(). */
    void setRcp(const RCP<PointerType>& ptr) 
    {ptr_=rcp_const_cast<const PointerType>(ptr);}
    /** \brief Protected non-const access to the underlying smart pointer. 
     *
     * This will be called by the nonConstPtr() method of the non-const Handle
     * subclass */
    RCP<PointerType> nonConstPtr() const
    {return rcp_const_cast<PointerType>(ptr_);}
  private:
    /** \brief . */
    RCP<const PointerType> ptr_;
  };

  /** \brief Generic templated handle class.
   *
   * In writing derived types, it is usually simplest to use the
   * TEUCHOS_HANDLE_CTORS macro to generate boilerplate constructor code.
   * 
   * There are two modes of construction:
   * construction from an existing RCP,
   * \code
   * RCP<Base> r = rcp(new Derived(blahblah));
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
    /** \brief . */
    Handle()
      : ConstHandle<PointerType>() {}
    /** \brief Construct with an existing RCP */
    Handle(const RCP<PointerType>& smartPtr)
      : ConstHandle<PointerType>() 
    {
      /* \brief We need to set the rcp in the base class */
      setRcp(smartPtr);
    }
    /** \brief Construct with a raw pointer to a Handleable.
     *
     * This will make a call to rcp() internally, thus removing that call from the
     * user interface.
     */
    explicit Handle(Handleable<PointerType>* rawPtr)
      : ConstHandle<PointerType>()
    {
      /* \brief We need to set the rcp in the base class. */
      setRcp(rawPtr->getRcp());
    }
    /** \brief Read/write access to the underlying smart pointer.
     */
    RCP<PointerType> ptr() const {return this->nonConstPtr();}
    /** \brief Access to non-const raw pointer. */
    PointerType* rawPtr() const {return this->nonConstPtr().get();}
  };

} // namespace Teuchos

/** \brief This helper macro defines boilerplate constructors for classes
 * deriving from Handle.
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
handle() : Teuchos::Handle<contents >() {;} \
handle(Teuchos::Handleable<contents >* rawPtr) : Teuchos::Handle<contents >(rawPtr) {;} \
handle(const Teuchos::RCP<contents >& smartPtr) : Teuchos::Handle<contents >(smartPtr){;}

/** \brief. This helper macro defines boilerplate constructors for classes
 * deriving from ConstHandle.
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
handle( Teuchos::ENull _null = Teuchos::null ) : Teuchos::ConstHandle<contents >() {;} \
handle(const Teuchos::ConstHandleable<contents >* rawPtr) : Teuchos::ConstHandle<contents >(rawPtr) {;} \
handle(const Teuchos::RCP<const contents >& smartPtr) : Teuchos::ConstHandle<contents >(smartPtr){;}

#endif // TEUCHOS_CONSTHANDLE_HPP
