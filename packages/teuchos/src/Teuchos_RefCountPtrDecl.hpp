// ///////////////////////////////////////////////////////////////////////
// RefCountPtr_decl.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef REFCOUNTPTR_DECL_H
#define REFCOUNTPTR_DECL_H

#include "Teuchos_ConfigDefs.hpp"

#ifdef REFCOUNTPTR_INLINE_FUNCS
#define REFCOUNTPTR_INLINE inline
#else
#define REFCOUNTPTR_INLINE
#endif

namespace Teuchos {

namespace PrivateUtilityPack {
	class RefCountPtr_node;
}

/** @defgroup RefCountPtr_stuff Reference counting smart pointer class for  automatic garbage collection.
 * \ingroup Misc_grp
 *
 * ToDo: Put in link to latex paper and shorten this discussion (perhaps just the quickstart).
 *
 * The testing program <tt>\ref TestRefCountPtr_grp "TestRefCountPtrMain.cpp"</tt> runs
 * <tt>RefCountPtr</tt> through the paces and shows some how it is used.
 */
//@{

/// Used to initialize a \c RefCountPtr object to NULL using an implicit conversion!
enum ENull { null };

///
/** Function class for deallocator that uses \c delete to delete a pointer which is used by \c RefCountPtr.
 */
template<class T>
class DeallocDelete
{
public:
	/// Deallocates a pointer \c ptr using \c delete ptr.
	void free( T* ptr ) { if(ptr) delete ptr; }
}; // end class DeallocDelete

///
/** Templated class for reference counted smart pointers.
 *
 * This is a class for smart reference counted pointer objects
 * that deletes an object (if the object is owned) that was allocated
 * by <tt>new</tt> after all refereces to it have been removed.
 *
 * For casting, non-member template functions are \ref RefCountPtr_conv "supplied"
 * to allow the equaivalient to implicit pointer casts and const casts.  Dynamic
 * casting with multiple inheritance is also handled by this implementation.
 *
 * The client must be very careful how this class is used.
 * Conversions from <tt>RefCountPtr<T1></tt> to
 * <tt>RefCountPtr<T2></tt> objects must be handled by the conversion
 * functions.  Never use <tt>get()</tt> to perform a conversion.
 *
 * For a more detailed discussion see this \ref RefCountPtr_stuff "description".
 */
template<class T>
class RefCountPtr {
public:
	///
	typedef T	element_type;
	///
	/** Initialize to NULL.
	 *
	 * This allows clients to write code like:
	 \code
	 RefCountPtr<int> p = null;
	 \endcode
	 or
	 \code
	 RefCountPtr<int> p;
	 \endcode
	 * and construct to <tt>NULL</tt>
	 */
	RefCountPtr( ENull null_arg = null );
	///
	/** Initialize from another <tt>RefCountPtr<T></tt> object.
	 *
	 * After construction, <tt>this</tt> and <tt>r_ptr</tt> will
	 * reference the same object.
	 *
	 * This form of the copy constructor is required even though the
	 * below more general templated version is sufficient since some
	 * compilers will generate this function automatically which will
	 * give an incorrect implementation.
	 *
	 * Postconditons:<ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	RefCountPtr(const RefCountPtr<T>& r_ptr);
	///
	/** Initialize from another <tt>RefCountPtr<T2></tt> object (implicit conversion only).
	 *
	 * This function allows the implicit conversion of smart pointer objects just
	 * like with raw C++ pointers.  Note that this function will only compile
	 * if the statement <tt>T1 *ptr = r_ptr.get()</tt> will compile.
	 *
	 * Postconditons: <ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	template<class T2>
	RefCountPtr(const RefCountPtr<T2>& r_ptr);
	///
	/** Removes a reference to a dynamically allocated object and possibly deletes
	 * the object if owned.
	 *
	 * Peforms <tt>delete ...</tt> if <tt>this->has_ownership() == true</tt> and
	 * <tt>this->count() == 1</tt>.  If <tt>this->count() == 1</tt> but <tt>this->has_ownership() == false</tt>
	 * then the object is not deleted.
	 * If <tt>this->count() > 1</tt> then then internal reference count
	 * shared by all the other related <tt>RefCountPtr<...></tt> objects for this shared
	 * object is deincremented by one.  If <tt>this->get() == NULL</tt> then nothing happens.
	 */
	~RefCountPtr();
	///
	/** Copy the pointer to the referenced object and increment the reference count.
	 *
	 * If <tt>this->has_ownership() == true</tt> and <tt>this->count() == 1</tt> before this operation
	 * is called, then the object pointed to by <tt>this->get()</tt> will be deleted (using <tt>delete</tt>)
	 * prior to binding to the pointer (possibly <tt>NULL</tt>) pointed to in <tt>r_ptr</tt>.
	 * Assignment to self (i.e. <tt>this->get() == r_ptr.get()</tt>) is harmless and this
	 * function does nothing.
	 *
	 * Postconditons:
	 * <ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	RefCountPtr<T>& operator=(const RefCountPtr<T>& r_ptr);
	///
	/** Pointer (<tt>-></tt>) access to members of underlying object.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>this->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
	 * </ul>
	 */
	T* operator->() const;
	///
	/** Dereference the underlying object.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>this->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
	 * </ul>
	 */
	T& operator*() const;
	///
    /** Get the raw C++ pointer to the underlying object.
	 */
	T* get() const;
	///
	/** Release the ownership of the underlying dynamically allocated object.
	 *
	 * After this function is called then the client is responsible for calling
	 * delete on the returned pointer no matter how many <tt>ref_count_prt<T></tt> objects
	 * have a reference to it.  If <tt>this-></tt>get() <tt>== NULL</tt>, then this call is
	 * meaningless.
	 *
	 * Note that this function does not have the exact same semantics as does
	 * <tt>auto_ptr<T>::release()</tt>.  In <tt>auto_ptr<T>::release()</tt>, <tt>this</tt>
	 * is set to <tt>NULL</tt> while here in RefCountPtr<T>:: release() only an ownership flag is set
	 * and <tt>this</tt> still points to the same object.  It would be difficult to duplicate
	 * the behavior of <tt>auto_ptr<T>::release()</tt> for this class.
	 *
	 * Postconditions:
	 * <ul>
	 * <li> <tt>this->has_ownership() == false</tt>
	 * </ul>
	 *
	 * @return Returns the value of <tt>this->get()</tt>
	 */
	T* release();
	///
	/** Return the number of <tt>RefCountPtr<></tt> objects that have a reference
	 * to the underlying pointer that is being shared.
	 *
	 * @return  If <tt>this->get() == NULL</tt> then this function returns 0.
	 * Otherwise, this function returns <tt>> 0</tt>.
	 */
	int count() const;
	///
	/** Call to give <tt>this</tt> and other <tt>RefCountPtr<></tt> objects ownership 
	 * of the referenced object <tt>this->get()</tt>.
	 *
	 * See ~RefCountPtr() above.  This function
	 * does nothing if <tt>this->get() == NULL</tt>.
	 *
	 * Postconditions:
	 * <ul>
	 * <li> If <tt>this->get() == NULL</tt> then
	 *   <ul>
	 *   <li> <tt>this->has_ownership() == false</tt> (always!).
	 *   </ul>
	 * <li> else
	 *   <ul>
	 *   <li> <tt>this->has_ownership() == true</tt>
	 *   </ul>
	 * </ul>
	 */
	void set_has_ownership();
	///
	/** Returns true if <tt>this</tt> has ownership of object pointed to by <tt>this->get()</tt> in order to delete it.
	 *
	 * See ~RefCountPtr() above.
	 *
	 * @return If this->get() <tt>== NULL</tt> then this function always returns <tt>false</tt>.
	 * Otherwise the value returned from this function depends on which function was
	 * called most recently, if any; set_has_ownership() (<tt>true</tt>)
	 * or release() (<tt>false</tt>).
	 */
	bool has_ownership() const;
	///
	/** Return if the the smart pointers share the same resource.
	 *
	 * This method does more than just check if <tt>this->get() == r_ptr.get()</tt>.
	 * It also checks to see if the underlying reference counting machinary is the
	 * same.
	 */
	bool shares_resource(const RefCountPtr<T>& r_ptr) const;

private:

	// //////////////////////////////////////
	// Private types

	typedef PrivateUtilityPack::RefCountPtr_node			node_t;

	// //////////////////////////////////////////////////////////////
	// Private data members

	T       *ptr_;  // NULL if this pointer is null
	node_t	*node_;	// NULL if this pointer is null

	// /////////////////////////////////////
	// Private member functions

	void assert_not_null() const;

public:
#ifndef DOXYGEN_COMPILE
	// These constructors should be private but I have not had good luck making
	// this portable (i.e. using friendship etc.) in the past
	RefCountPtr( T* p, bool has_ownership );
	template<class Dealloc_T>
	RefCountPtr( T* p, Dealloc_T dealloc, bool has_ownership );
	// This is a very bad breach of encapsulation that is needed since MS VC++ 5.0 will
	// not allow me to declare template functions as friends.
	RefCountPtr( T* p, node_t* node);
	T*&           access_ptr();
	node_t*&      access_node();
	node_t*       access_node() const;
#endif

};	// end class RefCountPtr<...>

/** Return a \c RefCountPtr object properly typed.
 *
 * The client should only create a RefCountPtr object given a pointer
 * return by \c new by calling this function.
 *
 * KL 10/07/03: The original code sensibly made owns_mem an optional
 * argument that defaulted to true, but that did not compile on Solaris.
 * I've changed owns_mem to a mandatory argument and created a second
 * method with no owns_mem argument to handle the default case.
 */
template<class T>
RefCountPtr<T> rcp( T* p );

///
/** Return a \c RefCountPtr object properly typed.
 *
 * The client should only create a RefCountPtr object given a pointer
 * return by \c new by calling this function.
 *
 * KL 10/07/03: The original code sensibly made owns_mem an optional
 * argument that defaulted to true, but that did not compile on Solaris.
 * I've changed owns_mem to a mandatory argument and created a second
 * method with no owns_mem argument to handle the default case.
 */
template<class T>
RefCountPtr<T> rcp( T* p, bool owns_mem);

///
/** Initialize from a raw pointer with a deallocation policy.
 *
 * @param  p       [in] Raw C++ pointer that \c this will represent.
 * @param  dealloc [in] Deallocator policy object (copied by value) that defines
 *                 a function <tt>void Dealloc_T::free(T* p)</tt> that will
 *                 free the underlying object.
 * @param  has_ownership
 *                 [in] If true then <tt>return</tt> is allowed to delete
 *                 the underlying pointer by calling <tt>dealloc.free(p)</tt>.
 *                 when all references have been removed.
 *
 * Postconditions:<ul>
 * <li> The function <tt>void Dealloc_T::free(T* p)</tt> exists.
 * </ul>
 *
 * Postconditions:<ul>
 * <li> <tt>return.get() == p</tt>
 * <li> If <tt>p == NULL</tt> then
 *   <ul>
 *   <li> <tt>return.count() == 0</tt>
 *   <li> <tt>return.has_ownership() == false</tt>
 *   </ul>
 * <li> else
 *   <ul>
 *   <li> <tt>return.count() == 1</tt>
 *   <li> <tt>return.has_ownership() == owns_mem</tt>
 *   </ul>
 * </ul>
 *
 * By default, <tt>return</tt> has ownership to delete the object
 * pointed to by <tt>p</tt> when <tt>return</tt> is deleted (see
 * <tt>~RefCountPtr())</tt>.  It is vitually important that if
 * <tt>owns_mem == true</tt> that the address <tt>p</tt> that is
 * passed in is the same address that was returned by <tt>new</tt>.
 * With multiple inheritance this is not always the case.  See the
 * above discussion.  This class is templated to accept a deallocator
 * object that will free the pointer.  The other functions use a
 * default deallocator of type <tt>DeallocDelete</tt> which has a method
 * <tt>DeallocDelete::free()</tt> which just calls <tt>delete p</tt>.
 */
template<class T, class Dealloc_T>
RefCountPtr<T> rcp( T* p, Dealloc_T dealloc, bool owns_mem );

///
/** Implicit cast from T1* to T2*.
 *
 * The function will compile only if (<tt>T2* p2 = p1.get();</tt>) compiles.
 *
 * This is to be used for conversions up an inheritance hierarchy and from non-const to
 * const and any other standard implicit pointer conversions allowed by C++.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_implicit_cast(const RefCountPtr<T1>& p1);

///
/** static cast from T1* to T2*.
 *
 * The function will compile only if (<tt>static_cast<T2*>(p1.get());</tt>) compiles.
 *
 * This can safely be used for conversion down an inheritance hierarchy
 * with polymorphic types only if <tt>dynamic_cast<T2>(p1.get()) == static_cast<T2>(p1.get())</tt>.
 * If not then you have to use #rcp_dynamic_cast<tt><T2>(p1)</tt>.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_static_cast(const RefCountPtr<T1>& p1);

///
/** Constant cast from T1* to T2*.
 *
 * This function will compile only if (<tt>const_cast<T2*>(p1.get());</tt>) compiles.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_const_cast(const RefCountPtr<T1>& p1);

///
/** Dynamic cast from T1* to T2*.
 *
 * This function will compile only if (<tt>dynamic_cast<T2*>(p1.get());</tt>) compiles.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_dynamic_cast(const RefCountPtr<T1>& p1);

//@}

/** \defgroup TestRefCountPtr_grp Testing program for RefCountPtr<>.
 *
 * \include TestRefCountPtrMain.cpp
 */
//@{
//@}

} // end namespace Teuchos

#endif	// REFCOUNTPTR_DECL_H
