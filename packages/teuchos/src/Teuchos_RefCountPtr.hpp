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

#ifndef TEUCHOS_REFCOUNTPTR_H
#define TEUCHOS_REFCOUNTPTR_H

/*! \file Teuchos_RefCountPtr.hpp
    \brief Reference-counted pointer class and non-member templated function implementations.
*/

#include "Teuchos_RefCountPtrDecl.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_map.hpp"

// /////////////////////////////////////////////////////////////////////////
// Inline implementations below, not for the client to look at.

namespace Teuchos {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace PrivateUtilityPack {

// Assert that the pointer is not null
void throw_null( const std::string &type_name );

// Node class to keep track of the delete address and
// the reference count for RefCountPtr<...>
class RefCountPtr_node {
public:
	RefCountPtr_node(bool has_ownership)
		: count_(1), has_ownership_(has_ownership), extra_data_map_(NULL)
	{}
	virtual ~RefCountPtr_node()
	{
		if(extra_data_map_) delete extra_data_map_;
	}
	int count() const {
		return count_;	
	}
	int incr_count() {
		return ++count_;
	}
	int deincr_count() {
		return --count_;
	}
	void has_ownership(bool has_ownership) {
		has_ownership_ = has_ownership;
	}
	bool has_ownership() const {
		return has_ownership_;
	}
	void set_extra_data( const any &extra_data, const std::string& name, bool force_unique );
	any& get_extra_data( const std::string& type_name, const std::string& name );
	const any& get_extra_data( const std::string& type_name, const std::string& name ) const {
		return const_cast<RefCountPtr_node*>(this)->get_extra_data(type_name,name);
	}
private:
	typedef Teuchos::map<std::string,any>  extra_data_map_t;
	int                 count_;
	bool                has_ownership_;
	extra_data_map_t    *extra_data_map_;
	// This is made a pointer to reduce overhead for the general case
	// where this is not used
	RefCountPtr_node();
	RefCountPtr_node(const RefCountPtr_node&);
	RefCountPtr_node& operator=(const RefCountPtr_node&);
};	// end class RefCountPtr_node;

// Implementation class for actually deleting the object if has_ownership() == true.
template<class T, class Dealloc_T>
class RefCountPtr_node_tmpl : public RefCountPtr_node {
public:

	//
	RefCountPtr_node_tmpl(T* p, Dealloc_T dealloc, bool has_ownership)
		: RefCountPtr_node(has_ownership), ptr_(p), dealloc_(dealloc)
	{}
	//
	Dealloc_T& get_dealloc() { return dealloc_; }
	//
	const Dealloc_T& get_dealloc() const { return dealloc_; }
	//
	~RefCountPtr_node_tmpl() {
		if( has_ownership() )
			dealloc_.free(ptr_);
	}

private:

	T           *ptr_;
	Dealloc_T   dealloc_;
	// not defined and not to be called
	RefCountPtr_node_tmpl();
	RefCountPtr_node_tmpl(const RefCountPtr_node_tmpl&);
	RefCountPtr_node_tmpl& operator=(const RefCountPtr_node_tmpl&);

}; // end class RefCountPtr_node_tmpl<T>

}	// end namespace PrivateUtilityPack 

#endif

// /////////////////////////////////////////////////////////////////////////////////
// Inline member functions for RefCountPtr<...>.

template<class T>
inline
RefCountPtr<T>::RefCountPtr( ENull )
	: ptr_(NULL)
	, node_(NULL)
{}

template<class T>
REFCOUNTPTR_INLINE
RefCountPtr<T>::RefCountPtr(const RefCountPtr<T>& r_ptr)
	: ptr_(r_ptr.ptr_), node_(r_ptr.node_)
{
	if(node_) node_->incr_count();
}

template<class T>
REFCOUNTPTR_INLINE
template <class T2>
RefCountPtr<T>::RefCountPtr(const RefCountPtr<T2>& r_ptr)
	: ptr_(const_cast<T2*>(r_ptr.get()))                 // will not compile if T1 is not an ancestor of T2
	, node_(const_cast<node_t*>(r_ptr.access_node()))
{
	if(node_) node_->incr_count();
}

template<class T>
REFCOUNTPTR_INLINE
RefCountPtr<T>::~RefCountPtr()
{
	if(node_ && node_->deincr_count() == 0 ) delete node_;
}

template<class T>
REFCOUNTPTR_INLINE
RefCountPtr<T>& RefCountPtr<T>::operator=(const RefCountPtr<T>& r_ptr) {
	if(node_) {
		if( r_ptr.node_ == node_ )
			return *this; // Assignment to self!
		if( !node_->deincr_count() ) {
			delete node_;
		}
	}
	ptr_  = r_ptr.ptr_;
	node_ = r_ptr.node_;
	if(node_) node_->incr_count();
	return *this;
}

template<class T>
inline
T* RefCountPtr<T>::operator->() const {
	assert_not_null();
	return ptr_;
}

template<class T>
inline
T& RefCountPtr<T>::operator*() const {
	assert_not_null();
	return *ptr_;
}

template<class T>
inline
T* RefCountPtr<T>::get() const {
	return ptr_;
}

template<class T>
REFCOUNTPTR_INLINE
T* RefCountPtr<T>::release() {
	if(node_)
		node_->has_ownership(false);
	return ptr_;
}

template<class T>
REFCOUNTPTR_INLINE
int RefCountPtr<T>::count() const {
	if(node_)
		return node_->count();
	return 0;
}

template<class T>
REFCOUNTPTR_INLINE
void RefCountPtr<T>::set_has_ownership() {
	if(node_)
		node_->has_ownership(true);
}

template<class T>
REFCOUNTPTR_INLINE
bool RefCountPtr<T>::has_ownership() const {
	if(node_)
		return node_->has_ownership();
	return false;
}

template<class T>
REFCOUNTPTR_INLINE
bool RefCountPtr<T>::shares_resource(const RefCountPtr<T>& r_ptr) const {
	return node_ == r_ptr.node_;
}

// private

template<class T>
inline
void RefCountPtr<T>::assert_not_null() const {
	if(!ptr_) PrivateUtilityPack::throw_null(typeid(T).name());
}

// very bad public functions

template<class T>
inline
RefCountPtr<T>::RefCountPtr( T* p, bool has_ownership )
	: ptr_(p)
	, node_( p ? new PrivateUtilityPack::RefCountPtr_node_tmpl<T,DeallocDelete<T> >(p,DeallocDelete<T>(),has_ownership) : NULL )
{}

template<class T>
REFCOUNTPTR_INLINE
template<class Dealloc_T>
RefCountPtr<T>::RefCountPtr( T* p, Dealloc_T dealloc, bool has_ownership )
	: ptr_(p)
	, node_( p ? new PrivateUtilityPack::RefCountPtr_node_tmpl<T,Dealloc_T>(p,dealloc,has_ownership) : NULL )
{}

template<class T>
inline
RefCountPtr<T>::RefCountPtr( T* p, node_t* node)
	: ptr_(p), node_(node)
{
	if(node_) node_->incr_count();
}

template<class T>
inline
T*& RefCountPtr<T>::access_ptr()
{	return ptr_; }

template<class T>
inline
typename RefCountPtr<T>::node_t*& RefCountPtr<T>::access_node()
{	return node_; }

template<class T>
inline
typename RefCountPtr<T>::node_t* RefCountPtr<T>::access_node() const
{	return node_; }

}	// end namespace Teuchos

// /////////////////////////////////////////////////////////////////////////////////
// Inline member functions for conversions for RefCountPtr<...>.

template<class T>
inline
Teuchos::RefCountPtr<T>
Teuchos::rcp( T* p, bool owns_mem )
{
	return RefCountPtr<T>(p,owns_mem);
}

template<class T, class Dealloc_T>
inline
Teuchos::RefCountPtr<T>
Teuchos::rcp( T* p, Dealloc_T dealloc, bool owns_mem )
{
	return RefCountPtr<T>(p,dealloc,owns_mem);
}

template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RefCountPtr<T2>
Teuchos::rcp_implicit_cast(const RefCountPtr<T1>& p1)
{
	T2 *check = p1.get();	// Make the compiler check if the conversion is legal
	RefCountPtr<T2> p2;
	if(p1.access_node()) {
		p2.access_ptr()  = check;
		p2.access_node() = const_cast<RefCountPtr<T1>&>(p1).access_node();
		p2.access_node()->incr_count();
	}
	return p2;
}

template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RefCountPtr<T2>
Teuchos::rcp_static_cast(const RefCountPtr<T1>& p1)
{
	T2 *check = static_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
	RefCountPtr<T2> p2;
	if(p1.access_node()) {
		p2.access_ptr()  = check;
		p2.access_node() = const_cast<RefCountPtr<T1>&>(p1).access_node();
		p2.access_node()->incr_count();
	}
	return p2;
}

template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RefCountPtr<T2>
Teuchos::rcp_const_cast(const RefCountPtr<T1>& p1)
{
	T2 *check = const_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
	RefCountPtr<T2> p2;
	if(p1.access_node()) {
		p2.access_ptr()  = check;
		p2.access_node() = const_cast<RefCountPtr<T1>&>(p1).access_node();
		p2.access_node()->incr_count();
	}
	return p2;
}

template<class T2, class T1>
REFCOUNTPTR_INLINE
Teuchos::RefCountPtr<T2>
Teuchos::rcp_dynamic_cast(const RefCountPtr<T1>& p1)
{
	RefCountPtr<T2> p2; // NULL by default
	if( p1.get() ) {
		T2 *check = dynamic_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
		if(check) {
			p2.access_ptr()  = check;
			p2.access_node() = const_cast<RefCountPtr<T1>&>(p1).access_node();
			p2.access_node()->incr_count();
		}
	}
	return p2;
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
void Teuchos::set_extra_data( const T1 &extra_data, const std::string& name, Teuchos::RefCountPtr<T2> *p, bool force_unique )
{
	*(*p); // Assert not NULL
	p->access_node()->set_extra_data( extra_data, name, force_unique );
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
T1& Teuchos::get_extra_data( RefCountPtr<T2>& p, const std::string& name )
{
	*p; // Assert not NULL
	return any_cast<T1>(p.access_node()->get_extra_data(typeid(T1).name(),name));
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1& Teuchos::get_extra_data( const RefCountPtr<T2>& p, const std::string& name )
{
	*p; // Assert not NULL
	return any_cast<T1>(p.access_node()->get_extra_data(typeid(T1).name(),name));
}

template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
Dealloc_T&
Teuchos::get_dealloc( RefCountPtr<T>& p )
{
	*p; // Assert not NULL
	PrivateUtilityPack::RefCountPtr_node_tmpl<typename Dealloc_T::ptr_t,Dealloc_T>
		*dnode = dynamic_cast<PrivateUtilityPack::RefCountPtr_node_tmpl<typename Dealloc_T::ptr_t,Dealloc_T>*>(p.access_node());
	TEST_FOR_EXCEPTION(
		dnode==NULL, std::logic_error
		,"get_dealloc<" << typeid(Dealloc_T).name() << "," << typeid(T).name() << ">(p): "
		<< "Error, requested type \'" << typeid(PrivateUtilityPack::RefCountPtr_node_tmpl<typename Dealloc_T::ptr_t,Dealloc_T>).name()
		<< "\' does not match actual type of the node \'" << typeid(*p.access_node()).name() << "!"
		);
	return dnode->get_dealloc();
}

template<class Dealloc_T, class T>
inline
const Dealloc_T& 
Teuchos::get_dealloc( const Teuchos::RefCountPtr<T>& p )
{
	return get_dealloc<Dealloc_T>(const_cast<RefCountPtr<T>&>(p));
}

#endif	// REFCOUNTPTR_H
