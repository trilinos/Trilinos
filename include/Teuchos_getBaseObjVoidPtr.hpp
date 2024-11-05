// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_GET_BASE_OBJ_VOID_PTR_HPP
#define TEUCHOS_GET_BASE_OBJ_VOID_PTR_HPP


#include "TeuchosCore_ConfigDefs.hpp"
#include "Teuchos_ConfigDefs.hpp"


#if defined(HAVE_TEUCHOSCORE_CXX11)
#  define HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR 1
#  include <type_traits>
#elif defined(HAVE_TEUCHOSCORE_BOOST) && defined(HAVE_TEUCHOSCORE_BOOST_IS_POLYMORPHIC)
#  define HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR 1
#  include <boost/type_traits/is_polymorphic.hpp>
#endif


namespace Teuchos {


#ifdef HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR


template<bool isPolymorphic, typename T>
class GetBaseObjVoidPtrImpl {};

template<typename T>
class GetBaseObjVoidPtrImpl<true, T> {
public:
  static const void* getBaseObjVoidPtr(T *p)
    {
      return dynamic_cast<const void*>(p);
    }
};


template<typename T>
class GetBaseObjVoidPtrImpl<false, T> {
public:
  static const void* getBaseObjVoidPtr(T *p)
    {
      return static_cast<const void*>(p);
    }
};


/** \brief Return a const void* pointing to the base of an object.
 *
 * This function uses the boost::is_polymorphic traits class to determine if
 * type T is a polymorphic type or a non-polymorphic and then calls
 * dynamic_cast or static_cast, respectively, to return the base pointer to
 * the object.
 *
 * The base pointer to an object is crtical to know if you need to determine
 * if two pointers are pointing to the same object or not.
 *
 * NOTE: This function will not even be defined if
 * HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR is not defined (which currently requires
 * boost support but that may change later).
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
const void* getBaseObjVoidPtr(T *p)
{
#if defined(HAVE_TEUCHOSCORE_CXX11)
  const bool isPoly = std::is_polymorphic<T>::value;
#elif defined(HAVE_TEUCHOSCORE_BOOST) && defined(HAVE_TEUCHOSCORE_BOOST_IS_POLYMORPHIC)
  const bool isPoly = boost::is_polymorphic<T>::value;
#endif
  typedef GetBaseObjVoidPtrImpl<isPoly, T> GBOVPT;
  return GBOVPT::getBaseObjVoidPtr(p);
}


#endif // HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR


}	// end namespace Teuchos


#endif // TEUCHOS_GET_BASE_OBJ_VOID_PTR_HPP
