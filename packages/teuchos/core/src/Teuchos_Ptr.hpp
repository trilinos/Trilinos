// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PTR_HPP
#define TEUCHOS_PTR_HPP


#include "Teuchos_PtrDecl.hpp"
#include "Teuchos_RCP.hpp"


namespace Teuchos {


namespace PtrPrivateUtilityPack {
TEUCHOSCORE_LIB_DLL_EXPORT void throw_null( const std::string &type_name );
} // namespace PtrPrivateUtilityPack


template<class T> inline
Ptr<T>::Ptr( ENull /*null_in*/ )
  : ptr_(0)
{}


template<class T> inline
Ptr<T>::Ptr( T *ptr_in )
  :ptr_(ptr_in)
{}


template<class T> inline
Ptr<T>::Ptr(const Ptr<T>& ptr_in)
  :ptr_(ptr_in.ptr_)
#ifdef TEUCHOS_DEBUG
  ,rcp_(ptr_in.access_rcp())
#endif
{}


template<class T>
template<class T2> inline
Ptr<T>::Ptr(const Ptr<T2>& ptr_in)
  :ptr_(ptr_in.get())
#ifdef TEUCHOS_DEBUG
  ,rcp_(ptr_in.access_rcp())
#endif
{}


template<class T> inline
Ptr<T>& Ptr<T>::operator=(const Ptr<T>& ptr_in)
{
  ptr_ = ptr_in.get();
#ifdef TEUCHOS_DEBUG
  rcp_ = ptr_in.rcp_;
#endif
  return *this;
}


template<class T> inline
T* Ptr<T>::operator->() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T> inline
T& Ptr<T>::operator*() const
{
  debug_assert_not_null();
  debug_assert_valid_ptr();
  return *ptr_;
}


template<class T> inline
T* Ptr<T>::get() const
{
  debug_assert_valid_ptr();
  return ptr_;
}


template<class T> inline
T* Ptr<T>::getRawPtr() const
{
  return get();
}


template<class T> inline
const Ptr<T>& Ptr<T>::assert_not_null() const
{
  if(!ptr_)
    PtrPrivateUtilityPack::throw_null(TypeNameTraits<T>::name());
  return *this;
}


template<class T> inline
bool Ptr<T>::is_null () const {
  return ptr_ == NULL;
}


template<class T> inline
const Ptr<T> Ptr<T>::ptr() const
{
  return *this;
}


template<class T> inline
Ptr<const T> Ptr<T>::getConst() const
{
  return ptr_implicit_cast<const T>(*this);
}


template<class T> inline
void Ptr<T>::debug_assert_valid_ptr() const
{
#ifdef TEUCHOS_DEBUG
  rcp_.access_private_node().assert_valid_ptr(*this);
#endif
}


#ifdef TEUCHOS_DEBUG


template<class T> inline
Ptr<T>::Ptr( const RCP<T> &p )
  : ptr_(p.getRawPtr()), rcp_(p)
{}


#endif // TEUCHOS_DEBUG


} // namespace Teuchos


template<class T>
std::ostream& Teuchos::operator<<( std::ostream& out, const Ptr<T>& p )
{
  out
    << TypeNameTraits<RCP<T> >::name() << "{"
    << "ptr="<<(const void*)(p.get()) // I can't find any alternative to this C cast :-(
    <<"}";
  return out;
}


#endif // TEUCHOS_PTR_HPP
