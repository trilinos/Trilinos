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

#ifndef TEUCHOS_REFCOUNTPTR_DECL_H
#define TEUCHOS_REFCOUNTPTR_DECL_H

/*! \file Teuchos_RefCountPtrDecl.hpp
    \brief Reference-counted pointer class and non-member templated function
	definitions
*/

#include "Teuchos_any.hpp"

#ifdef REFCOUNTPTR_INLINE_FUNCS
#define REFCOUNTPTR_INLINE inline
#else
#define REFCOUNTPTR_INLINE
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace MemMngPack {} // ToDo: Take out latter!
#endif

/** \class Teuchos::DeallocDelete
    \brief Policy class for deallocator that uses <tt>delete</tt> to delete a
    pointer which is used by <tt>RefCountPtr</tt>.
 */

/** \class Teuchos::RefCountPtr
    \brief Templated class for reference counted smart pointers.
 *
 * This is a class for smart reference counted pointer objects
 * that deletes an object (if the object is owned) that was allocated
 * by <tt>new</tt> after all refereces to it have been removed.
 *
 * For casting, non-member template functions are \ref RefCountPtr_stuff "supplied"
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

namespace Teuchos {

namespace PrivateUtilityPack {
	class RefCountPtr_node;
}

/** \defgroup RefCountPtr_stuff Reference counting smart pointer class for automatic garbage collection.
 * 
 * ToDo: Put in link to latex paper and shorten this discussion (perhaps just the quickstart).
 *
 */
//@{

/// Used to initialize a <tt>RefCountPtr</tt> object to NULL using an implicit conversion!
enum ENull { null };

///
template<class T>
class DeallocDelete
{
public:
	/// Gives the type (required)
	typedef T ptr_t;
	/// Deallocates a pointer <tt>ptr</tt> using <tt>delete ptr</tt> (required).
	void free( T* ptr ) { if(ptr) delete ptr; }
};

///
template<class T>
class RefCountPtr {
public:
	///
	typedef T	element_type;
	///
	/** \brief Initialize <tt>RefCountPtr<T></tt> to NULL.
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
	/** \brief Initialize from another <tt>RefCountPtr<T></tt> object.
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
	/** \brief Initialize from another <tt>RefCountPtr<T2></tt> object (implicit conversion only).
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
	/** \brief Removes a reference to a dynamically allocated object and possibly deletes
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
	/** \brief Copy the pointer to the referenced object and increment the reference count.
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
	/** \brief Pointer (<tt>-></tt>) access to members of underlying object.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>this->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
	 * </ul>
	 */
	T* operator->() const;
	///
	/** \brief Dereference the underlying object.
	 *
	 * Preconditions:
	 * <ul>
	 * <li> <tt>this->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
	 * </ul>
	 */
	T& operator*() const;
	///
    /** \brief Get the raw C++ pointer to the underlying object.
	 */
	T* get() const;
	///
	/** \brief Release the ownership of the underlying dynamically allocated object.
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
	/** \brief Return the number of <tt>RefCountPtr<></tt> objects that have a reference
	 * to the underlying pointer that is being shared.
	 *
	 * @return  If <tt>this->get() == NULL</tt> then this function returns 0.
	 * Otherwise, this function returns <tt>> 0</tt>.
	 */
	int count() const;
	///
	/** \brief Give <tt>this</tt> and other <tt>RefCountPtr<></tt> objects ownership 
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
	/** \brief Returns true if <tt>this</tt> has ownership of object pointed to by <tt>this->get()</tt> in order to delete it.
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
	/** \brief Returns true if the the smart pointers share the same resource.
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

///
/** \brief Create a <tt>RefCountPtr</tt> object properly typed.
 *
 * @param  p  [in] Pointer to an object to be reference counted.
 * @param owns_mem
 *            [in] If <tt>owns_mem==true</tt>  then <tt>delete p</tt>
 *            will be called when the last reference to this object
 *            is removed.  If <tt>owns_mem==false</tt> then nothing
 *            will happen to delete the the object pointed to by
 *            <tt>p</tt> when the last reference is removed.
 *
 * Preconditions:<ul>
 * <li> If <tt>owns_mem==true</tt> then <tt>p</tt> must have been
 *      created by calling <tt>new</tt> to create the object since
 *      <tt>delete p</tt> will be called eventually.
 * </ul>
 *
 * If the pointer <tt>p</tt> did not come from <tt>new</tt> then
 * either the client should use the version of <tt>rcp()</tt> that
 * that uses a deallocator policy object or should pass in 
 * <tt>owns_mem = false</tt>.
 */
template<class T>
RefCountPtr<T> rcp( T* p, bool owns_mem
#ifndef __sun
	= true
#endif
	);
#ifdef __sun // RAB: 20040303: Sun needs to fixe there &^**%F$ compiler
template<class T> inline RefCountPtr<T> rcp( T* p ) { return rcp(p,true); }
#endif

///
/** \brief Initialize from a raw pointer with a deallocation policy.
 *
 * @param  p       [in] Raw C++ pointer that \c this will represent.
 * @param  dealloc [in] Deallocator policy object (copied by value) that defines
 *                 a function <tt>void Dealloc_T::free(T* p)</tt> that will
 *                 free the underlying object.
 * @param  owns_mem
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
/** \brief Implicit cast of underlying <tt>RefCountPtr</tt> type from T1* to T2*.
 *
 * The function will compile only if (<tt>T2* p2 = p1.get();</tt>) compiles.
 *
 * This is to be used for conversions up an inheritance hierarchy and from non-const to
 * const and any other standard implicit pointer conversions allowed by C++.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_implicit_cast(const RefCountPtr<T1>& p1);

///
/** \brief Static cast of underlying <tt>RefCountPtr</tt> type from T1* to T2*.
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
/** \brief Constant cast of underlying <tt>RefCountPtr</tt> type from T1* to T2*.
 *
 * This function will compile only if (<tt>const_cast<T2*>(p1.get());</tt>) compiles.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_const_cast(const RefCountPtr<T1>& p1);

///
/** \brief Dynamic cast of underlying <tt>RefCountPtr</tt> type from T1* to T2*.
 *
 * This function will compile only if (<tt>dynamic_cast<T2*>(p1.get());</tt>) compiles.
 */
template<class T2, class T1>
RefCountPtr<T2> rcp_dynamic_cast(const RefCountPtr<T1>& p1);

///
/** \brief Set extra data associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  extra_data
 *               [in] Data object that will be set (copied)
 * @param  name  [in] The name given to the extra data.  The value of
 *               <tt>name</tt> together with the data type <tt>T1</tt> of the
 *               extra data must be unique from any other such data or
 *               the other data will be overwritten.
 * @param  p     [out] On output, will be updated with the input <tt>extra_data</tt>
 * @param  force_unique
 *               [in] Determines if this type and name pair must be unique
 *               in which case if an object with this same type and name
 *               already exists, then an exception will be thrown.
 *               The default is <tt>true</tt> for safety.
 *
 * If there is a call to this function with the same type of extra
 * data <tt>T1</tt> and same arguments <tt>p</tt> and <tt>name</tt>
 * has already been made, then the current piece of extra data already
 * set will be overwritten with <tt>extra_data</tt>.  However, if the
 * type of the extra data <tt>T1</tt> is different, then the extra
 * data can be added and not overwrite existing extra data.  This
 * means that extra data is keyed on both the type and name.  This
 * helps to minimize the chance that clients will unexpectedly
 * overwrite data by accident.
 *
 * When the last <tt>RefcountPtr</tt> object is removed the underlying
 * reference-counted object is deleted before any of the extra data
 * that has been associated with this object.  The extra data objects
 * will then be destoried in a first-in basis.  In other words, the
 * first extra data object added will be deleted first, the second
 * extra data object will be deleted second and so on.  This must be
 * considered when multiple pieces of extra data are being added if
 * the order of distruction is significant.
 *
 * Preconditions:<ul>
 * <li> <tt>p->get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> If this function has already been called with the same template
 *      type <tt>T1</tt> for <tt>extra_data</tt> and the same string <tt>name</tt>
 *      and <tt>force_unique==true</tt>, then an <tt>std::invalid_argument</tt>
 *      exception will be thrown.
 * </ul>
 *
 * Note, this function is made a non-member function to be consistent
 * with the non-member <tt>get_extra_data()</tt> functions.
 */
template<class T1, class T2>
void set_extra_data( const T1 &extra_data, const std::string& name, RefCountPtr<T2> *p, bool force_unique
#ifndef __sun
	 = true
#endif
	);
#ifdef __sun
template<class T1, class T2>
inline void set_extra_data( const T1 &extra_data, const std::string& name, RefCountPtr<T2> *p )
{ set_extra_data( extra_data, name, p, true ); }
#endif

///
/** \brief Get a non-const reference to extra data associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  p    [in] Smart pointer object that extra data is being extraced from.
 * @param  name [in] Name of the extra data.
 *
 * @return Returns a non-const reference to the extra_data object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> <tt>name</tt> and <tt>T1</tt> must have been used in a previous
 *      call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>).
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 */
template<class T1, class T2>
T1& get_extra_data( RefCountPtr<T2>& p, const std::string& name );

///
/** \brief Get a const reference to extra data associated with a <tt>RefCountPtr</tt> object.
 *
 * @param  p    [in] Smart pointer object that extra data is being extraced from.
 * @param  name [in] Name of the extra data.
 *
 * @return Returns a const reference to the extra_data object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> <tt>name</tt> and <tt>T1</tt> must have been used in a previous
 *      call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>).
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * Also note that this const version is a false sense of security
 * since a client can always copy a const <tt>RefCountPtr</tt> object
 * into a non-const object and then use the non-const version to
 * change the data.  However, its presence will help to avoid some
 * types of accidental changes to this extra data.
 */
template<class T1, class T2>
const T1& get_extra_data( const RefCountPtr<T2>& p, const std::string& name );

///
/** \brief Return a non-<tt>const</tt> reference to the underlying deallocator object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> The deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 *      (throws <tt>std::logic_error</tt>)
 * </ul>
 *
 */
template<class Dealloc_T, class T>
Dealloc_T& get_dealloc( RefCountPtr<T>& p );

///
/** \brief Return a <tt>const</tt> reference to the underlying deallocator object.
 *
 * Preconditions:<ul>
 * <li> <tt>p.get() != NULL</tt> (throws <tt>std::logic_error</tt>)
 * <li> The deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 *      (throws <tt>std::logic_error</tt>)
 * </ul>
 *
 * Note that the <tt>const</tt> version of this function provides only
 * a very ineffective attempt to avoid accidental changes to the
 * deallocation object.  A client can always just create a new
 * non-<tt>const</tt> <tt>RefCountPtr<T></tt> object from any
 * <tt>const</tt> <tt>RefCountPtr<T></tt> object and then call the
 * non-<tt>const</tt> version of this function.
 */
template<class Dealloc_T, class T>
const Dealloc_T& get_dealloc( const RefCountPtr<T>& p );

//@}

} // end namespace Teuchos

#endif	// TEUCHOS_REFCOUNTPTR_DECL_H
