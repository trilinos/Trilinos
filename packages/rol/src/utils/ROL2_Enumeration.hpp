// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL2_ENUMERATION_HPP
#define ROL2_ENUMERATION_HPP

#include "ROL2_Types.hpp"

/** \file  ROL2_Enumeration.hpp
 *  \brief Defines a standardized ROL2 scoped enum type and utility functions
 *
 */

namespace ROL2 {

constexpr bool is_enum_v = std::is_enum<T>::value;

template<class ENumType, class ReturnType>
using enable_if_enum_t = std::enable_if_t<is_enum_v<ENumType>,ReturnType>;

/** 
 *  \brief Universal valid value function for all ROL2 enums
 *         All enums in ROL2 follow the same design
 *
 *          enum class ENumType : SomeIntegralType {
 *            FirstValue = 0,
 *            ...
 *            Last 
 *          };
 *   
 *         The ENumType values always start at 0, are consecutive,
 *         and the largest value is denoted as Last. Values are 
 *         considered valid if they are at least 0 and less than Last.
 *         This function will throw an exception if the value is less than
 *         0 or greater than ENumType::Last
 */

template<class ENumType>
enable_if_enum_t<ENumType,bool>
is_valid_enum_value( const ENumType& e ) {
  using integral_type = std::underlying_type_t<ENumType>;
  constexpr auto first   = static_cast<integral_type>(0);
  constexpr auto last    = static_cast<integral_type>(ENumType::Last);
            auto current = static_cast<integral_type>(e);
  ROL_TEST_FOR_EXCEPTION( (current < first) || (current > last),
                          std::out_of_range, "The enum value is out of range" );
  return (first <= current) && (current < last);
};

/** \brief Universal post-increment operator of ROL2 class enums
*/
template<typename ENumType>
enable_if_enum_t<ENumType,ENumType&>
operator++ ( ENumType& e ) {
  using int_type = std::underlying_type_t<ENumType>;
  auto e_int = static_cast<int_type>(e);
  constexpr auto last = static_cast<int_type>(ENumType::Last);
  ROL_TEST_FOR_EXCEPTION( e_int >= last, std::out_of_range, 
                          "The enum value is out of range" );
  e = static_cast<ENumType>(++e_int);
  return e; 
}

/** \brief Universal pre-increment operator of ROL2 class enums
*/
template<typename ENumType>
enable_if_enum_t<ENumType,ENumType&>
operator++ ( ENumType& e, int ) {
  using int_type = std::underlying_type_t<ENumType>;
  auto e_int = static_cast<int_type>(e);
  constexpr auto last = static_cast<int_type>(ENumType::Last);
  ROL_TEST_FOR_EXCEPTION( e_int >= last, std::out_of_range, 
                          "The enum value is out of range" );
  e = static_cast<ENumType>(++e_int);
  return e; 
}

template<class ENumType>
enable_if_enum_t<ENumType,std::ostream&>
operator << ( std::ostream& os, ENumType e ) {
  using int_type = std::underlying_type_t<ENumType>;
  os << static_cast<int_type>(e);
  return os;
}

template<class ENumType, class RHS>
enable_if_enum_t<ENumType,bool>
operator == ( ENumType e, RHS rhs ) {
  using int_type = std::underlying_type_t<ENumType>;
  return static_cast<int_type>(e) == rhs;
}

template<class LHS, class ENumType>
enable_if_enum_t<ENumType,bool>
operator == ( LHS lhs, ENumType e ) {
  using int_type = std::underlying_type_t<ENumType>;
  return static_cast<int_type>(e) == lhs;
}

template<class ENumType, class RHS>
enable_if_enum_t<ENumType,bool>
operator > ( ENumType e, RHS rhs ) {
  using int_type = std::underlying_type_t<ENumType>;
  return static_cast<int_type>(e) > rhs;
}

template<class LHS, class ENumType>
enable_if_enum_t<ENumType,bool>
operator > ( LHS lhs, ENumType e ) {
  using int_type = std::underlying_type_t<ENumType>;
  return lhs > static_cast<int_type>(e);
}

template<class ENumType, class RHS>
enable_if_enum_t<ENumType,bool>
operator < ( ENumType e, RHS rhs ) {
  using int_type = std::underlying_type_t<ENumType>;
  return static_cast<int_type>(e) < rhs;
}

template<class LHS, class ENumType>
enable_if_enum_t<ENumType,bool>
operator < ( LHS lhs, ENumType e ) {
  using int_type = std::underlying_type_t<ENumType>;
  return lhs < static_cast<int_type>(e);
}

template<class ENumType, class RHS>
enable_if_enum_t<ENumType,bool>
operator >= ( ENumType e, RHS rhs ) {
  using int_type = std::underlying_type_t<ENumType>;
  return static_cast<int_type>(e) >= rhs;
}

template<class LHS, class ENumType>
enable_if_enum_t<ENumType,bool>
operator >= ( LHS lhs, ENumType e ) {
  using int_type = std::underlying_type_t<ENumType>;
  return lhs >= static_cast<int_type>(e);
}

template<class ENumType, class RHS>
enable_if_enum_t<ENumType,bool>
operator <= ( ENumType e, RHS rhs ) {
  using int_type = std::underlying_type_t<ENumType>;
  return static_cast<int_type>(e) <= rhs;
}

template<class LHS, class ENumType>
enable_if_enum_t<ENumType,bool>
operator <= ( LHS lhs, ENumType e ) {
  using int_type = std::underlying_type_t<ENumType>;
  return lhs <= static_cast<int_type>(e);
}



} // namespace ROL2

#endif // ROL2_ENUMERATION_HPP

