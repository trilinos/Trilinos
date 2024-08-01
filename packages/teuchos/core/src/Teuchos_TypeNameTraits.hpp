// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_TYPE_NAME_TRAITS_HPP_
#define _TEUCHOS_TYPE_NAME_TRAITS_HPP_

/*! \file Teuchos_TypeNameTraits.hpp
 \brief Defines basic traits returning the
    name of a type in a portable and readable way.
*/

// mfh 30 Jan 2013: Thanks to Jim Willenbring for reporting this, and
// to Mike Glass and Paul Lin for updating the fix for dealing with a
// bug in IBM's XL C++ compiler.  The update was necessary due to a
// relapse of the bug in a newer version of the compiler.
//
// If you don't have this update, you can fix the problem by defining
// the macro TEUCHOS_TYPE_NAME_TRAITS_OLD_IBM when compiling anything
// that includes this header file.  If you have the current version of
// this file, then you don't need to do anything.
#if defined(__IBMCPP__) && ( __IBMCPP__ < 900 || __IBMCPP__ == 1210 )
# define TEUCHOS_TYPE_NAME_TRAITS_OLD_IBM
#endif

#include <typeinfo>

#include "Teuchos_ConfigDefs.hpp"

namespace  Teuchos {


/** \brief Demangle a C++ name if valid.
 *
 * The name must have come from <tt>typeid(...).name()</tt> in order to be
 * valid name to pass to this function.
 *
 * \ingroup teuchos_language_support_grp
 */
TEUCHOSCORE_LIB_DLL_EXPORT std::string demangleName( const std::string &mangledName );


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
#ifndef TEUCHOS_TYPE_NAME_TRAITS_OLD_IBM
  static std::string concreteName( const T& t )
#else
  // the IBM compilers on AIX have a problem with const
  static std::string concreteName( T t )
#endif
    {
      return demangleName(typeid(t).name());
    }
};


/** \brief Template function for returning the concrete type name of a
 * passed-in object.
 *
 * Uses the traits class TypeNameTraits so the behavior of this function can
 * be specialized in every possible way.  The default return value is
 * typically derived from <tt>typeid(t).name()</tt>.
 *
 * \ingroup teuchos_language_support_grp
 */
template<typename T>
std::string typeName( const T &t )
{
  typedef typename std::remove_const_t<T> ncT;
#ifndef TEUCHOS_TYPE_NAME_TRAITS_OLD_IBM
  return TypeNameTraits<ncT>::concreteName(t);
#else
  // You can't pass general objects to AIX by value as above.  This means that
  // you will not get the concrete name printed on AIX but that is life on
  // such compilers.
  return TypeNameTraits<ncT>::name();
#endif
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
  typedef typename std::remove_const_t<T> ncT;
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
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(signed char);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(unsigned char);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(short int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(long int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(unsigned short int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(unsigned int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(unsigned long int);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(float);
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(double);

#ifdef HAVE_TEUCHOSCORE_QUADMATH
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(__float128);
#endif // HAVE_TEUCHOSCORE_QUADMATH

#ifdef HAVE_TEUCHOS_LONG_DOUBLE
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(long double);
#endif // HAVE_TEUCHOS_LONG_DOUBLE

template<typename T>
class TEUCHOSCORE_LIB_DLL_EXPORT TypeNameTraits<T*> {
public:
  typedef T* T_ptr;
  static std::string name() { return TypeNameTraits<T>::name() + "*"; }
  static std::string concreteName(T_ptr) { return name(); }
};


template<>
class TEUCHOSCORE_LIB_DLL_EXPORT TypeNameTraits<std::string> {
public:
  static std::string name() { return "string"; }
  static std::string concreteName(const std::string&)
    { return name(); }
};


template<>
class TEUCHOSCORE_LIB_DLL_EXPORT TypeNameTraits<void*> {
public:
  static std::string name() { return "void*"; }
  static std::string concreteName(const std::string&) { return name(); }
};

// mfh 31 Jul 2012: Specialization for "void" will hopefully fix
// compile errors on Windows, such as the following:
//
// http://testing.sandia.gov/cdash/viewBuildError.php?buildid=611137
//
// I'm imitating the specialization of void* above.
template<>
class TEUCHOSCORE_LIB_DLL_EXPORT TypeNameTraits<void> {
public:
  static std::string name() { return "void"; }
  static std::string concreteName(const std::string&) { return name(); }
};


#ifdef HAVE_TEUCHOS_COMPLEX


template<typename T>
class TEUCHOSCORE_LIB_DLL_EXPORT TypeNameTraits<std::complex<T> > {
public:
  static std::string name()
    { return "complex<"+TypeNameTraits<T>::name()+">"; }
  static std::string concreteName(const std::complex<T>&)
    { return name(); }
};


#endif // HAVE_TEUCHOS_COMPLEX



} // namespace Teuchos


#endif // _TEUCHOS_TYPE_NAME_TRAITS_HPP_
