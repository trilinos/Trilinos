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

#ifndef TEUCHOS_REFCOUNTPTR_HPP
#define TEUCHOS_REFCOUNTPTR_HPP

/*! \file Teuchos_RefCountPtr.hpp
    \brief Reference-counted pointer class and non-member templated function implementations.
*/
/** \example example/RefCountPtr/cxx_main.cpp
    This is an example of how to use the <tt>Teuchos::RefCountPtr</tt> class.
*/
/** \example test/RefCountPtr/cxx_main.cpp
    This is a more detailed testing program that uses all of the <tt>Teuchos::RefCountPtr</tt> class.
*/

#include "Teuchos_RefCountPtrDecl.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_map.hpp"

//#define TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES // Define this on command line to keep track of this!

// /////////////////////////////////////////////////////////////////////////
// Inline implementations below, not for the client to look at.

namespace Teuchos {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace PrivateUtilityPack {

//
class RefCountPtr_node;

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
	void set_extra_data( const any &extra_data, const std::string& name, EPrePostDestruction destroy_when, bool force_unique );
	any& get_extra_data( const std::string& type_name, const std::string& name );
	const any& get_extra_data( const std::string& type_name, const std::string& name ) const {
		return const_cast<RefCountPtr_node*>(this)->get_extra_data(type_name,name);
	}
	any* get_optional_extra_data( const std::string& type_name, const std::string& name );
	const any* get_optional_extra_data( const std::string& type_name, const std::string& name ) const {
		return const_cast<RefCountPtr_node*>(this)->get_optional_extra_data(type_name,name);
	}
protected:
  void pre_delete_extra_data() {
    if(extra_data_map_) impl_pre_delete_extra_data();
  }
private:
  struct extra_data_entry_t {
    extra_data_entry_t() : destroy_when(POST_DESTROY) {}
    extra_data_entry_t( const any &_extra_data, EPrePostDestruction _destroy_when )
      : extra_data(_extra_data), destroy_when(_destroy_when) {}
    any extra_data;
    EPrePostDestruction destroy_when;
  };  
	typedef Teuchos::map<std::string,extra_data_entry_t> extra_data_map_t;
	int                 count_;
	bool                has_ownership_;
	extra_data_map_t    *extra_data_map_;
	// Above is made a pointer to reduce overhead for the general case
	// where this is not used
  void impl_pre_delete_extra_data();
  // Not defined and not to be called
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
    this->pre_delete_extra_data();
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

// Add new RCP to global list
void add_new_RefCountPtr_node( RefCountPtr_node* rcp_node, const std::string &info );

// Remove RCP from global list
void remove_RefCountPtr_node( RefCountPtr_node* rcp_node );

// Print global list
void print_active_RefCountPtr_nodes(std::ostream &out);

// Print global list on destruction
class PrintActiveRefCountPtrNodes {
public:
  PrintActiveRefCountPtrNodes();
  ~PrintActiveRefCountPtrNodes();
  void foo();
private:
  static int count_;
};

}	// end namespace PrivateUtilityPack 

} // namespace Teuchos

#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES

namespace {
// This static variable should be delcared before all other static variables
// that depend on RefCountPtr and therefore This static varaible should be
// deleted *after* all of these other static variables that depend on
// RefCountPtr go away!
Teuchos::PrivateUtilityPack::PrintActiveRefCountPtrNodes printActiveRefCountPtrNodes;
} // namespace
#endif

namespace Teuchos {

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
	if(node_ && node_->deincr_count() == 0 ) {
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
    printActiveRefCountPtrNodes.foo(); // Make sure this object is used!
    remove_RefCountPtr_node(node_);
#endif
    delete node_;
  }
}

template<class T>
REFCOUNTPTR_INLINE
RefCountPtr<T>& RefCountPtr<T>::operator=(const RefCountPtr<T>& r_ptr)
{
  if( this == &r_ptr )
    return *this; // Assignment to self
	if( node_ && !node_->deincr_count() ) {
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
    remove_ArrayRefCountPtr_node(node_);
#endif
    delete node_;
  }
	ptr_  = r_ptr.ptr_;
	node_ = r_ptr.node_;
	if(node_) node_->incr_count();
	return *this;
}

template<class T>
inline
T* RefCountPtr<T>::operator->() const {
#ifdef TEUCHOS_REFCOUNTPTR_ASSERT_NONNULL
	assert_not_null();
#endif
	return ptr_;
}

template<class T>
inline
T& RefCountPtr<T>::operator*() const {
#ifdef TEUCHOS_REFCOUNTPTR_ASSERT_NONNULL
	assert_not_null();
#endif
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
template <class T2>
bool RefCountPtr<T>::shares_resource(const RefCountPtr<T2>& r_ptr) const {
	return node_ == r_ptr.access_node();
  // Note: above, r_ptr is *not* the same class type as *this so we can not
  // access its node_ member directly!  This is an interesting detail to the
  // C++ protected/private protection mechanism!
}

template<class T>
inline
const RefCountPtr<T>& RefCountPtr<T>::assert_not_null() const {
	if(!ptr_) PrivateUtilityPack::throw_null(typeid(T).name());
	return *this;
}

// very bad public functions

template<class T>
inline
RefCountPtr<T>::RefCountPtr( T* p, bool has_ownership )
	: ptr_(p)
	, node_( p ? new PrivateUtilityPack::RefCountPtr_node_tmpl<T,DeallocDelete<T> >(p,DeallocDelete<T>(),has_ownership) : NULL )
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
  if(node_) {
    std::ostringstream os;
    os << "{T=\'"<<typeid(T).name()<<"\',Concrete T=\'"<<typeid(*p).name()<<"\',p="<<p<<",has_ownership="<<has_ownership<<"}";
    add_new_RefCountPtr_node(node_,os.str());
  }
#endif
}

template<class T>
REFCOUNTPTR_INLINE
template<class Dealloc_T>
RefCountPtr<T>::RefCountPtr( T* p, Dealloc_T dealloc, bool has_ownership )
	: ptr_(p)
	, node_( p ? new PrivateUtilityPack::RefCountPtr_node_tmpl<T,Dealloc_T>(p,dealloc,has_ownership) : NULL )
{
#ifdef TEUCHOS_SHOW_ACTIVE_REFCOUNTPTR_NODES
  if(node_) {
    std::ostringstream os;
    os << "{T=\'"<<typeid(T).name()<<"\',Concrete T=\'"<<typeid(*p).name()<<"\',p="<<p<<",has_ownership="<<has_ownership<<"}";
    add_new_RefCountPtr_node(node_,os.str());
  }
#endif
}

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
// Inline non-member functions for RefCountPtr

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

template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::is_null( const RefCountPtr<T> &p )
{
  return p.get() == NULL;
}

template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const RefCountPtr<T> &p, ENull )
{
  return p.get() == NULL;
}

template<class T>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const RefCountPtr<T> &p, ENull )
{
  return p.get() != NULL;
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator==( const RefCountPtr<T1> &p1, const RefCountPtr<T2> &p2 )
{
  return p1.access_node() == p2.access_node();
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
bool Teuchos::operator!=( const RefCountPtr<T1> &p1, const RefCountPtr<T2> &p2 )
{
  return p1.access_node() != p2.access_node();
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
Teuchos::rcp_dynamic_cast(const RefCountPtr<T1>& p1, bool throw_on_fail)
{
	RefCountPtr<T2> p2; // NULL by default
	if( p1.get() ) {
		T2 *check = NULL;
		if(throw_on_fail)
			check = &dyn_cast<T2>(*p1);
		else
			check = dynamic_cast<T2*>(p1.get()); // Make the compiler check if the conversion is legal
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
void Teuchos::set_extra_data( const T1 &extra_data, const std::string& name, Teuchos::RefCountPtr<T2> *p, EPrePostDestruction destroy_when, bool force_unique )
{
	p->assert_not_null();
	p->access_node()->set_extra_data( any(extra_data), name, destroy_when, force_unique );
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
T1& Teuchos::get_extra_data( RefCountPtr<T2>& p, const std::string& name )
{
	p.assert_not_null();
	return any_cast<T1>(p.access_node()->get_extra_data(typeid(T1).name(),name));
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1& Teuchos::get_extra_data( const RefCountPtr<T2>& p, const std::string& name )
{
	p.assert_not_null();
	return any_cast<T1>(p.access_node()->get_extra_data(typeid(T1).name(),name));
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
T1* Teuchos::get_optional_extra_data( RefCountPtr<T2>& p, const std::string& name )
{
	p.assert_not_null();
  any *extra_data = p.access_node()->get_optional_extra_data(typeid(T1).name(),name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}

template<class T1, class T2>
REFCOUNTPTR_INLINE
const T1* Teuchos::get_optional_extra_data( const RefCountPtr<T2>& p, const std::string& name )
{
	p.assert_not_null();
  any *extra_data = p.access_node()->get_optional_extra_data(typeid(T1).name(),name);
  if( extra_data ) return &any_cast<T1>(*extra_data);
  return NULL;
}

template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
Dealloc_T&
Teuchos::get_dealloc( RefCountPtr<T>& p )
{
	p.assert_not_null();
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

template<class Dealloc_T, class T>
REFCOUNTPTR_INLINE
Dealloc_T*
Teuchos::get_optional_dealloc( RefCountPtr<T>& p )
{
	p.assert_not_null();
	PrivateUtilityPack::RefCountPtr_node_tmpl<typename Dealloc_T::ptr_t,Dealloc_T>
		*dnode = dynamic_cast<PrivateUtilityPack::RefCountPtr_node_tmpl<typename Dealloc_T::ptr_t,Dealloc_T>*>(p.access_node());
  if(dnode)
    return &dnode->get_dealloc();
  return NULL;
}

template<class Dealloc_T, class T>
inline
const Dealloc_T*
Teuchos::get_optional_dealloc( const Teuchos::RefCountPtr<T>& p )
{
	return get_optional_dealloc<Dealloc_T>(const_cast<RefCountPtr<T>&>(p));
}

template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const RefCountPtr<T>& p )
{
  out
    << TypeNameTraits<RefCountPtr<T> >::name() << "{"
    << "ptr="<<(const void*)(p.get()) // I can't find any alternative to this C cast :-(
    <<",node="<<p.access_node()
    <<",count="<<p.count()
    <<"}";
  return out;
}

#endif // TEUCHOS_REFCOUNTPTR_HPP
