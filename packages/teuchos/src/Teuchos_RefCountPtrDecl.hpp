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

// Uncomment this if not supported or define it in arguments
// to the C++ compiler (even better)
#ifndef REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
#define REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
#endif

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
 * Described here is a template class (Teuchos::RefCountPtr<tt><...></tt>) and a set of template
 * functions for implementing automatic garbage collection in C++ using smart reference counted pointers.
 * This design is based partly on the interface for <tt>std::auto_ptr<T></tt> and
 * Items 28 and 29 in "More Effective C++" by Scott Myers.  In short,
 * using this class allows one client to dynamically create an object (using new),
 * pass it around to bunch of other clients which also need it and
 * then never requiring any client to explicitly call delete but the object
 * will still (almost magically) be deleted when everyone is finished using it.
 * Boy, doesn't this sound like garbage collection like is implimented
 * in Perl and Java.  However, to realize this wonderful behavior requires
 * following some rules.  These rules will be spelled out later.
 *
 * The following class hierarchy will be used to demonstrate this smart pointer design.
 \code

 class A { public: virtual ~A() {} A& operator=(const A&) {} virtual void f() {} };
 class B1 : virtual public A {};
 class B2 : virtual public A {};
 class C : virtual public B1, virtual public B2 {};
 class D {};
 class E : public D {};
 \endcode
 * Note that the classes <tt>A</tt>, <tt>B1</tt>, <tt>B2</tt> and <tt>C</tt> are polymorphic (with
 * multiple inheritance with virtual base classes to boot) while the classes <tt>D</tt> and <tt>E</tt> are not.  
 *
 * In the following description, all of the code examples are written as though there was a 
 * <tt>using namespace %Teuchos;</tt> declaration is the current scope.
 *
 * A smart reference counted pointer to a dynammicallly allocated object of type <tt>A</tt> is initialized like:
 \code

 RefCountPtr<A>             a_ptr   = rcp(new A);    // i.e. A       *       a_ptr   = new A;
 RefCountPtr<const A>       ca_ptr  = rcp(new A);    // i.e. const A *       ca_ptr  = new A;
 const RefCountPtr<A>       a_cptr  = rcp(new A);    // i.e. A       * const a_cptr  = new A;
 const RefCountPtr<const A> ca_cptr = rcp(new A);    // i.e. const A * const ca_cptr = new A;
 \endcode
 * The above code shows how all the various combinations for const or non-const pointers
 * to const or non-const objects are expressed.  There is no automatic conversion from
 * raw pointers to smart pointers for many good reasons.  To allow such an implicit
 * conversion almost guarentees memory allocation/deallocation errors will occur.
 * Therefore, the templated function \c rcp<T>(T* p) must be used to initialize a smart
 * pointer (with one exception explained below for initializing to \c NULL).
 * Through the magic of the member functions \c operator->() and \c operator*(), these referenced
 * objects can be accessed exactly the same as for raw pointers like:
 \code

 a_ptr->f();
 (*a_ptr).f();
 \endcode
 * Therefore, using a smart reference counted pointer is very similar to
 * using a raw C++ pointer.  Some things you can do with smart reference countered
 * pointers that are similar as for raw pointers are:
 \code

 RefCountPtr<A>
   a_ptr1     = rcp(new A),   // Initialize them from a raw pointer from new
   a_ptr2     = rcp(new A);   // ""
 A *ra_ptr1   = new A,        // ""
   *ra_ptr2   = new A;        // ""
 a_ptr1       = rcp(ra_ptr2); // Assign from a raw pointer (only do this once!)
 a_ptr2       = a_ptr1;       // Assign one smart pointer to another
 a_ptr1       = rcp(ra_ptr1); // Assign from a raw pointer (only do this once!)
 a_ptr1->f();                 // Access a member using ->
 ra_ptr1->f();                // ""
 *a_ptr1      = *a_ptr2;      // Dereference the objects and assign
 *ra_ptr1     = *ra_ptr2;     // "" 
 \endcode
 * What makes smart pointers different however is that the above piece of code
 * does not create any memory leaks that would have otherwise been the case
 * if <tt>a_ptr1</tt> and <tt>a_ptr2</tt> where raw C++ pointers.
 * However, these smart reference counted pointers can not be used everywhere a
 * raw pointer can be.  For instance the following statements will not compile:
 \code

 a_ptr1++;           // Error, pointer arithmetic ++, --, +, - etc. not defined!
 a_ptr1 == ra_ptr1;  // Error, comparision operators ==, !=, <=, >= etc. not defined!
 \endcode
 * Since a smart reference counted pointer can not be used everywhere a raw C++ pointer
 * can be, there is a means for getting at the raw C++ pointer with the <tt>get()</tt> function
 * (same as with \c std::auto_ptr<T>::get()).  For example, to check if two smart pointer
 * objects contain references to the same object you would check:
 \code

 if( a_ptr1.get() == a_ptr2.get() )
   std::cout << "a_ptr1 and a_ptr2 point to the same object\n";
 \endcode
 * 
 * The conversion of smart pointers according to inheritance is just as easy
 * as with raw C++ pointers if member templates are supported.  This option is already
 * supported but must be enabled by defining the preprocessor macro
 * \c REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS.  For raw C++ pointers, the compiler will
 * implicitly cast up an inheritance hiearchy (i.e. <tt>*B1 => *A</tt> or <tt>*E => *D</tt>).
 * For smart reference counted pointers with member template functions supported, this is just
 * as easy.  For example, the below code compiles and runs just fine:
 \code

 RefCountPtr<C>  c_ptr  = rcp(new C);
 RefCountPtr<A>  a_ptr  = c_ptr;
 RefCountPtr<B1> b1_ptr = c_ptr;
 RefCountPtr<D>  d_ptr  = rcp(new E);
 \endcode
 * If member template functions are not supported (and the macro
 * \c REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS is not defined)
 * the conversion of smart pointers according to inheritance is not so
 * nice as with raw C++ pointers.  To accomplish casts that the C++
 * compiler would perform for you, you must use the following non-member template function:
 \code

 template<class T2, class T1>
 RefCountPtr<T2> rcp_implicit_cast(const RefCountPtr<T1>& p1);
 \endcode
 * The above code without member template functons would have to be written like:
 \code

 RefCountPtr<C>  c_ptr  = rcp(new C);
 RefCountPtr<A>  a_ptr  = rcp_implicit_cast<A>(c_ptr);
 RefCountPtr<B1> b1_ptr = rcp_implicit_cast<B1>(c_ptr);
 RefCountPtr<D>  d_ptr  = rcp_implicit_cast<D>(rcp(new E));
 \endcode
 * Okay, so having to use <tt>rcp_implicit_cast<...></tt> to explicitly perform all casts
 * that should be performed implicitly is a pain, but if your compiler does not support
 * member template functions then you have no choice.
 *
 * To perform other non-implicit type conversions with pointers such as
 * <tt>static_cast<...></tt>, <tt>const_cast<...></tt> and <tt>dynamic_cast<...></tt> there are
 * corresponding non-member template functions \c rcp_static_cast(), \c rcp_const_cast() and
 * \c rcp_dynamic_cast() respectively with the following prototypes:
 \code

 template<class T2, class T1>
 RefCountPtr<T2> rcp_static_cast(const RefCountPtr<T1>& p1);

 template<class T2, class T1>
 RefCountPtr<T2> rcp_const_cast(const RefCountPtr<T1>& p1);

 template<class T2, class T1>
 RefCountPtr<T2> rcp_dynamic_cast(const RefCountPtr<T1>& p1);
 \endcode
 * The usage of these conversion template functions looks very similar to the
 * syntax with raw C++ pointers.  For example:
 \code

 RefCountPtr<const C>   c_ptr   = rcp(new C);
 RefCountPtr<const A>   ca_ptr  = c_ptr;
 RefCountPtr<const C>   cc_ptr1 = rcp_dynamic_cast<const C>(ca_ptr);  // Safe and checked at runtime!
 RefCountPtr<const C>   cc_ptr2 = rcp_static_cast<const C>(ca_ptr);   // Potentially unsafe and unchecked!
 RefCountPtr<A>         a_ptr   = rcp_const_cast<A>(ca_ptr);          // Cast away const
 \endcode
 * Just as with the case with the built-in C++ conversion operators, some types
 * of conversions will not compile.  For example:
 \code

 RefCountPtr<C>          c_ptr1 = rcp(new C);
 RefCountPtr<A>          a_ptr  = c_ptr;
 RefCountPtr<const C>    c_ptr2 = rcp_dynamic_cast<const C>(a_ptr);  // Error, can't const and dyn cast at same time!
 RefCountPtr<const C>    c_ptr3 = rcp_dynamic_cast<const C>(
                                        const_cast<const A>(a_ptr)
                                      );                               // Okay
 RefCountPtr<D>          d_ptr  = rcp(new D);
 RefCountPtr<E>          e_ptr1 = rcp_dynamic_cast<E>(d_ptr);        // Error, D and E are not polymorphic!
 RefCountPtr<E>          e_ptr2 = rcp_static_cast<E>(d_ptr);         // Okay?
 \endcode
 * There are a few rules that you must absolutly follow when using
 * <tt>RefCountPtr<...></tt> objects:
 *
 * 1) Never give a raw pointer to more than one <tt>RefCountPtr<...></tt> object.
 * If you do you must call release on all but one of them or multiple calls
 * to delete will be performed and result in a runtime error.
 \code

 {
   C  *rc_ptr = new C;
   ...
   RefCountPtr<C> c_ptr1 = rcp(rc_ptr); // c_ptr1 has ownership to delete!
   ...
   RefCountPtr<C> c_ptr2 = rcp(rc_ptr); // c_ptr2 knows nothing of c_ptr1 but assumes has ownership to delete!
   ...
   // At end of block when c_ptr2 is destroyed, delete on rc_ptr will be called
   // then delete will be called again when c_ptr1 is destroyed!  This is
   // a runtime error that may crash the program.
 }
 \endcode
 *
 * 2) When sharing objects from polymorphic classes using multiple inheritance,
 * always make the first RefCountPtr<> object have the same type used with
 * \c new.
 \code

 {
   C*                rc_ptr1 = new C;
   RefCountPtr<B1> b1_ptr1(rc_ptr1); // Bad! (void*)b1_ptr1.get() != (void*)rc_ptr1 in general!
   C*                rc_ptr2 = new C;
   RefCountPtr<B1> b1_ptr2 = rcp(rc_ptr2); // Okay! rcp() will create a RefCountPtr<C> object
   // When delete is called on the internal shared pointer in the destructor for b1_ptr2,
   // the address (void*)rc_ptr2 will be used which is the same that was returned from new.
   // However, when delete is called on the internal shared pointer in the destructor
   // for b1_ptr1 the address (void*)(B1*)rc_ptr1 will be used which in general may not
   // be the same as (void*)rc_ptr2 which was returned from new.  Most C++ implementations
   // will (incorrectly) just call free() after the destructor is called which will yield
   // a runtime error since this not the same address returned from malloc().
 }
 \endcode
 * The above detail is very tricky but if you follow the guidelines below this will not
 * come up to bite you.
 *
 * 3) Always perform conversions using the non-member template functions described
 * above.  Never, and I mean never, perform a conversion using
 * \c Teuchos::RefCountPtr::get().
 \code

 {
   RefCountPtr<C> c_ptr = rcp(new C);       // c_ptr has ownership to detete
   ...
   RefCountPtr<A> a_ptr = rcp(c_ptr.get()); // This performs the conversion but
                                              // a_ptr knowns nothing of c_ptr.
    ...
	// At end of block when a_ptr is destroyed, delete will be called
	// then delete will be called again when c_ptr is destroyed!  This is
	// a runtime error that may crash the program.
 }
 \endcode
 * The above example is very similar to 1) in that two <tt>RefCountPtr<...></tt>
 * objects have ownership to delete a dynamically allocated object but do
 * know of the other <tt>RefCountPtr<...></tt> object.  Therefore when these
 * smart pointer objects are deleted they don't know that any other clients
 * still needed it so that it will not delete the object.
 *
 * 4) Never give any <tt>RefCountPtr<...></tt> objects the pointer to a static
 * (i.e. global) or automatic (i.e. created on the stack) object without
 * passing <tt>has_ownership==false</tt> to the constructor
 * Teuchos::RefCountPtr::RefCountPtr<tt>(p,has_ownership)</tt>
 * or calling calling Teuchos::RefCountPtr::release() immediatly after.
 * If you don't, then <tt>delete</tt> will be called on this pointer and cause a memory
 * coruption at runtime.
 \code

 {
   C c1; // Allocated by the compiler on the stack (not dynamically allocated!)
   C c2; // ""
   ...
   RefCountPtr<C> c1_ptr = &c1;      // Oh no! c1_ptr has ownership to delete!
   RefCountPtr<C> c2_ptr(&c2,false); // Okay! will not try to delete c2!
   ...
   // At the end of the block when c1_ptr is destroyed delete will be called
   // on &c1 even through it was not allocated with new.  This is a runtime
   // error that will probably crash the program.  However, when c2_ptr is destroyed
   // nothing (bad) will happen.
 }
 \endcode
 *
 * 5) Never give a raw C++ pointer returned from <tt>RefCountPtr<...>::get()</tt>
 * to an external entity which may delete it unless you call <tt>RefCountPtr<...>::release()</tt>
 * first.  It is probably better to have a memory leak than to perform multiple
 * deletes on the same pointer.
 \code

 void f( A* ra_ptr )
 {
   ...
   if( ... ) delete ra_ptr;
   ...
 }

 void main()
 {
   ...
   RefCountPtr<C> c_ptr = rcp(new C);
   ...
   f(c_ptr.get());  // May delete this pointer?
   ...
   // Reguardless of whether f(...) deleted the pointer returned by 
   // c_ptr.get() or not, c_ptr will call delete on this pointer when
   // it is destroyed at the end of the block.
 }
 \endcode
 *
 * To help not break the above hard and fast rules, the following guidelines
 * should be followed when using Teuchos::RefCountPtr<tt><...></tt>.
 *
 * 1) Get the raw pointer returned from <tt>new</tt> into a <tt>RefCountPtr<...></tt> object
 * as soon as possible.
 *
 * The simplest way to do this is to call <tt>new</tt> directly in the call the the 
 * nonmember function \c rcp() as shown below.
 \code
 RefCountPtr<BaseClass> ptr = rcp(new DerivedClass());
 \endcode
 * To initialize to NULL, use the implicit conversion from the dummy enum value \c null:
 \code
 RefCountPtr<BaseClass> ptr = null;
 \endcode
 * These are the only ways to initialize a \c RefCountPtr object with a raw
 * pointer without calling the explicit constructor and a casual user should never
 * call an explicit constructor directly.
 *
 * 2) Use <tt>RefCountPtr<...></tt> for all of your dynamic memory handling!  The more
 * consistently you use this class, the better your life will be.
 *
 * If everyone plays nice then <tt>RefCountPtr<...></tt> can give you the illusion of
 * garbage collection in C++.  However, it is not hard to screw things up if you
 * are not careful.
 *
 * The testing program <tt>\ref TestRefCountPtr_grp "TestRefCountPtrMain.cpp"</tt> runs
 * \c RefCountPtr through the paces and shows some how it is used.
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
 * The client must be very careful how this class is used.  Conversions
 * from <tt>RefCountPtr<T1></tt> to <tt>RefCountPtr<T2></tt> objects must be
 * handled by the conversion functions.  Do not use get() to perform
 * a conversion.
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
	 * This is equivalent to calling <tt>RefCountPtr(NULL)</tt> accept it
	 * allows clients to write code like:
	 \code
	 RefCountPtr<int> p = null;
	 \endcode
	 */
	RefCountPtr( ENull );

	///
	/** Equivalent to calling the constructor <tt>RefCountPtr(p,DeallocDelete(),has_ownership)</tt>.
	 */
	explicit RefCountPtr( T* p = NULL, bool has_ownership = true );
#ifdef REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
	///
	/** Initialize from a raw pointer (does not allow implicit conversions!).
	 *
	 * By default, <tt>this</tt> has ownership to delete the object pointed to
	 * by <tt>p</tt> when <tt>this</tt> is deleted (see ~RefCountPtr()).  It is
	 * vitually important that if <tt>has_ownership == true</tt> that the address
	 * <tt>p</tt> that is passed in is the same address that was returned by
	 * <tt>new</tt>.  With multiple inheritance this is not always the case.
	 * See the above discussion.  This constructor is templated to accept a
	 * deallocator object that will free pointer.  The other constructor
	 * uses yy default, the type \c DeallocDelete which has a method
	 * \c DeallocDelete::free() which just calls \c delete ptr.
	 *
	 * @param  p       [in] Raw C++ pointer that \c this will represent.
	 * @param  dealloc [in] Deallocator object (copied by value) that defines
	 *                 a function <tt>void Dealloc_T::free(T* p)</tt> that will
	 *                 free the underlying object.
	 * @param  has_ownership
	 *                 [in] If true then \c this will have be allowed to delete
	 *                 the underlying pointer by calling <tt>dealloc.free(p)</tt>.
	 *                 when all references have been removed.
	 *
	 * Postconditions:
	 * <ul>
	 * <li> The function <tt>void Dealloc_T::free(T* p)</tt> exists.
	 * <li> <tt>this->get() == p</tt>
	 * <li> If <tt>p == NULL</tt> then
	 *   <ul>
	 *   <li> <tt>this->count() == 0</tt>
	 *   <li> <tt>this->has_ownership() == false</tt>
	 *   </ul>
	 * <li> else
	 *   <ul>
	 *   <li> <tt>this->count() == 1</tt>
	 *   <li> <tt>this->has_ownership() == has_ownership</tt>
	 *   </ul>
	 * </ul>
	 *
	 * Warning, this function is declared \c explicit so don't expect an implicit
	 * conversion from a raw pointer to a \c RefCountPtr object.  Instead, use
	 * the function \c rcp() or the type \c null instead.
	 */
	template<class Dealloc_T>
	explicit RefCountPtr( T* p, Dealloc_T dealloc, bool has_ownership );
#endif
	///
	/** Initialize from another <tt>RefCountPtr<T></tt> object.
	 *
	 * After construction, <tt>this</tt> and <tt>r_ptr</tt> will reference the same
	 * object.
	 *
	 * Postconditons:
	 * <ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	RefCountPtr(const RefCountPtr<T>& r_ptr);
#ifdef REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
	///
	/** Initialize from another <tt>RefCountPtr<T2></tt> object (implicit conversion only).
	 *
	 * This function allows the implicit conversion of smart pointer objects just
	 * like with raw C++ pointers.  Note that this function will only compile
	 * if the statement <tt>T1 *ptr = r_ptr.get()</tt> will compile.
	 *
	 * Postconditons:
	 * <ul>
	 * <li> <tt>this->get() == r_ptr.get()</tt>
	 * <li> <tt>this->count() == r_ptr.count()</tt>
	 * <li> <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
	 * <li> If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
	 * </ul>
	 */
	template<class T2>
	RefCountPtr(const RefCountPtr<T2>& r_ptr);
#endif
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

	/// Assert that <tt>this->get() != NULL</tt>, throws <tt>std::logic_error</tt>.
	void assert_not_null() const;

private:

	// //////////////////////////////////////
	// Private types

	typedef PrivateUtilityPack::RefCountPtr_node			node_t;

	// //////////////////////////////////////////////////////////////
	// Private data members

	T       *ptr_;  // NULL if this pointer is null
	node_t	*node_;	// NULL if this pointer is null

public:
	// This is a very bad breach of encapsulation that is needed since MS VC++ 5.0 will
	// not allow me to declare template functions as friends.
#ifndef DOXYGEN_COMPILE
	RefCountPtr( T* p, node_t* node);
	T*&           access_ptr();
	node_t*&      access_node();
	node_t*       access_node() const;
#endif

};	// end class RefCountPtr<...>

///
/** Return a \c RefCountPtr object properly typed.
 *
 * The client should only create a RefCountPtr object given a pointer
 * return by \c new by calling this function.
 */
template<class T>
RefCountPtr<T> rcp( T* p, bool owns_mem = true );

#ifdef REFCOUNTPTR_TEMPLATE_CLASS_TEMPLATE_FUNCTIONS
///
/** Return a \c RefCountPtr object properly typed.
 *
 * The client should only create a RefCountPtr object given a pointer
 * return by \c new by calling this function.
 */
template<class T, class Dealloc_T>
RefCountPtr<T> rcp( T* p, Dealloc_T dealloc, bool owns_mem );
#endif

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
