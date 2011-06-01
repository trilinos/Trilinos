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

#ifndef TEUCHOS_GET_BASE_OBJ_VOID_PTR
#define TEUCHOS_GET_BASE_OBJ_VOID_PTR


#include "Teuchos_ConfigDefs.hpp"


#if defined(HAVE_TEUCHOS_BOOST) && defined(HAS_TEUCHOS_BOOST_IS_POLYMORPHIC)


// Internal Trilinos code should check if the macro is defined to see if
// getBaseObjPtr() is supported or not.
#define HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR 1


#include <boost/type_traits/is_polymorphic.hpp>


namespace Teuchos {


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
  typedef GetBaseObjVoidPtrImpl<boost::is_polymorphic<T>::value, T> GBOVPT;
  return GBOVPT::getBaseObjVoidPtr(p);
}


}	// end namespace Teuchos


#endif // defined(HAVE_TEUCHOS_BOOST) && defined(HAS_TEUCHOS_BOOST_IS_POLYMORPHIC)


#endif // TEUCHOS_GET_BASE_OBJ_VOID_PTR
