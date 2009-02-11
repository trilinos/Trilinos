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

#ifndef _TEUCHOS_TYPE_NAME_TRAITS_HPP_
#define _TEUCHOS_TYPE_NAME_TRAITS_HPP_

/*! \file Teuchos_TypeNameTraits.hpp
 \brief Defines basic traits returning the
    name of a type in a portable and readable way.
*/

#include "Teuchos_ConstTypeTraits.hpp"


namespace  Teuchos {


/** \brief Demangle a C++ name if valid.
 *
 * The name must have come from <tt>typeid(...).name()</tt> in order to be
 * valid name to pass to this function.
 *
 * \ingroup teuchos_language_support_grp
 */
std::string demangleName( const std::string &mangledName );


/** \brief Default traits class that just returns <tt>typeid(T).name()</tt>.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
class TypeNameTraits {
public:
  /** \brief . */
  static std::string name()
    {
      return demangleName(typeid(T).name());
    }
  /** \brief . */
#ifndef _AIX
  static std::string concreteName( const T& t )
#else
  // the IBM compilers on AIX have a problem with const
  static std::string concreteName( T& t )
#endif
    {
      return demangleName(typeid(t).name());
    }
};


/** \brief Template function for returning the type name of a passed-in
 * object.
 *
 * Uses the traits class TypeNameTraits so the behavior of this function can
 * be specialized in every possible way.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
std::string typeName( const T &t )
{
  typedef typename ConstTypeTraits<T>::NonConstType ncT;
  return TypeNameTraits<ncT>::name();
}


/** \brief Template function for returning the type name of the actual
 * concrete name of a passed-in object.
 *
 * Uses the traits class TypeNameTraits so the behavior of this function can
 * be specialized in every possible way.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
std::string concreteTypeName( const T &t )
{
  typedef typename ConstTypeTraits<T>::NonConstType ncT;
  return TypeNameTraits<ncT>::concreteName(t);
}


#define TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(TYPE) \
template<> \
class TypeNameTraits<TYPE> { \
public: \
  static std::string name() { return (#TYPE); } \
  static std::string concreteName(const TYPE&) { return name(); } \
} \

TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(bool);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(char);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(short int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(long int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(float);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(double);


template<typename T>
class TypeNameTraits<T*> {
public:
  typedef T* T_ptr;
  static std::string name() { return TypeNameTraits<T>::name() + "*"; }
  static std::string concreteName(T_ptr) { return name(); }
};


template<>
class TypeNameTraits<std::string> {
public:
  static std::string name() { return "string"; }
  static std::string concreteName(const std::string&)
    { return name(); }
};


template<>
class TypeNameTraits<void*> {
public:
  static std::string name() { return "void*"; }
  static std::string concreteName(const std::string&) { return name(); }
};


#ifdef HAVE_TEUCHOS_COMPLEX


template<typename T>
class TypeNameTraits<std::complex<T> > {
public:
  static std::string name()
    { return "complex<"+TypeNameTraits<T>::name()+">"; }
  static std::string concreteName(const std::complex<T>&)
    { return name(); }
};


#endif // HAVE_TEUCHOS_COMPLEX
 

} // namespace Teuchos


#endif // _TEUCHOS_TYPE_NAME_TRAITS_HPP_
