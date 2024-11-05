// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_ORDINALTRAITS_HPP_
#define _TEUCHOS_ORDINALTRAITS_HPP_

/*! \file Teuchos_OrdinalTraits.hpp
  \brief Defines basic traits for the ordinal field type
*/

#include "Teuchos_ConfigDefs.hpp"
#include <limits>

/*! \struct Teuchos::OrdinalTraits
  \brief This structure defines some basic traits for the ordinal field type.

  Ordinal traits are an essential part of templated codes.  This structure offers
  the basic traits of the templated ordinal type, like defining zero and one.

  For the general type, or default implementation, an aborting function
  is defined which should restrict implementations from using ordinal traits other than
  the defined specializations.

  \note The defined specializations for OrdinalTraits are: \c int and \c long \c int.
*/

/* This is the default structure used by OrdinalTraits<T> to produce a compile time
   error when the specialization does not exist for type <tt>T</tt>.
*/
namespace Teuchos {

template<class T>
struct UndefinedOrdinalTraits
{
  //! This function should not compile if there is an attempt to instantiate!
  static inline T notDefined() { return T::this_type_is_missing_a_specialization(); }
};

template<class T>
struct OrdinalTraits {

  //! Allows testing to see if ordinal traits machine parameters are defined.
  static const bool hasMachineParameters = false;

  //! Returns representation of zero for this ordinal type.
  static inline T zero()                     { return UndefinedOrdinalTraits<T>::notDefined(); }

  //! Returns representation of one for this ordinal type.
  static inline T one()                      { return UndefinedOrdinalTraits<T>::notDefined(); }

  //! \brief Returns a value designating the maximum value accessible by code using OrdinalTraits.
  /*! For a signed ordinal, this will typically be the maximum positive value. However, for an unsigned ordinal, this will typically be one less than the
   *  maximum value, because the maximum value will typically be reserved for designating an invalid value.
   */
  static inline T max()                  { return UndefinedOrdinalTraits<T>::notDefined(); }

  //! Returns a value designating an invalid number. For signed types, this is typically negative one; for unsigned types, it is typically the largest value.
  static inline T invalid()                  { return UndefinedOrdinalTraits<T>::notDefined(); }

  //! Returns name of this ordinal type.
  static inline std::string name()           { return UndefinedOrdinalTraits<T>::notDefined(); }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<>
struct OrdinalTraits<char> {
  static const bool hasMachineParameters = false;
  static inline char zero()                   {return(0);}
  static inline char one()                    {return(1);}
  static inline char invalid()                {return(std::numeric_limits<char>::max());}
  static inline char max()                    {return(std::numeric_limits<char>::max()-one());}
  static inline std::string name()            {return("char");}
};

template<>
struct OrdinalTraits<short int> {
  static const bool hasMachineParameters = false;
  static inline short int zero()              {return(0);}
  static inline short int one()               {return(1);}
  static inline short int invalid()           {return(-1);}
  static inline short int max()               {return(std::numeric_limits<short int>::max());}
  static inline std::string name()            {return("short int");}
};

template<>
struct OrdinalTraits<int> {
  static const bool hasMachineParameters = false;
  static inline int zero()                   {return(0);}
  static inline int one()                    {return(1);}
  static inline int invalid()                {return(-1);}
  static inline int max()                    {return(std::numeric_limits<int>::max());}
  static inline std::string name()           {return("int");}
};

template<>
struct OrdinalTraits<unsigned int> {
  static const bool hasMachineParameters = false;
  static inline unsigned int zero()                   {return(0);}
  static inline unsigned int one()                    {return(1);}
  static inline unsigned int invalid()                {return(std::numeric_limits<unsigned int>::max());}
  static inline unsigned int max()                    {return(std::numeric_limits<unsigned int>::max()-1);}
  static inline std::string name()           {return("unsigned int");}
};

template<>
struct OrdinalTraits<long int> {
  static const bool hasMachineParameters = false;
  static inline long int zero()              {return(static_cast<long int>(0));}
  static inline long int one()               {return(static_cast<long int>(1));}
  static inline long int invalid()           {return(static_cast<long int>(-1));}
  static inline long int max()               {return(std::numeric_limits<long int>::max());}
  static inline std::string name()           {return("long int");}
};

template<>
struct OrdinalTraits<long unsigned int> {
  static const bool hasMachineParameters = false;
  static inline long unsigned int zero()              {return(static_cast<long unsigned int>(0));}
  static inline long unsigned int one()               {return(static_cast<long unsigned int>(1));}
  static inline long unsigned int invalid()           {return(std::numeric_limits<long unsigned int>::max());}
  static inline long unsigned int max()               {return(std::numeric_limits<long unsigned int>::max()-1);}
  static inline std::string name()           {return("long unsigned int");}
};

template<>
struct OrdinalTraits<long long int> {
  static const bool hasMachineParameters = false;
  static inline long long int zero()          {return(static_cast<long long int>(0));}
  static inline long long int one()           {return(static_cast<long long int>(1));}
  static inline long long int invalid()       {return(static_cast<long long int>(-1));}
  static inline long long int max()           {return(std::numeric_limits<long long int>::max());}
  static inline std::string name()            {return("long long int");}
};

template<>
struct OrdinalTraits<unsigned long long int> {
  static const bool hasMachineParameters = false;
  static inline unsigned long long int zero()          {return(static_cast<unsigned long long int>(0));}
  static inline unsigned long long int one()           {return(static_cast<unsigned long long int>(1));}
  static inline unsigned long long int invalid()       {return(std::numeric_limits<unsigned long long int>::max());}
  static inline unsigned long long int max()           {return(std::numeric_limits<unsigned long long int>::max()-1);}
  static inline std::string name()                     {return("unsigned long long int");}
};

#ifdef HAVE_TEUCHOS___INT64

template<>
struct OrdinalTraits<__int64> {
  static const bool hasMachineParameters = false;
  static inline __int64 zero() {return(static_cast<__int64>(0));}
  static inline __int64 one() {return(static_cast<__int64>(1));}
  static inline __int64 invalid() {return(std::numeric_limits<__int64>::max());}
  static inline __int64 max() {return(std::numeric_limits<__int64>::max()-1);}
  static inline std::string name() {return("__int64");}
};

template<>
struct OrdinalTraits<unsigned __int64> {
  static const bool hasMachineParameters = false;
  static inline unsigned __int64 zero() {return(static_cast<unsigned __int64>(0));}
  static inline unsigned __int64 one() {return(static_cast<unsigned __int64>(1));}
  static inline unsigned __int64 invalid() {return(std::numeric_limits<unsigned __int64>::max());}
  static inline unsigned __int64 max() {return(std::numeric_limits<unsigned __int64>::max()-1);}
  static inline std::string name() {return("unsigned __int64");}
};

#endif // HAVE_TEUCHOS___INT64

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Teuchos

#endif // _TEUCHOS_ORDINALTRAITS_HPP_
