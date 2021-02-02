#pragma once
#ifndef ROL2_ENUMERATION_HPP
#define ROL2_ENUMERATION_HPP

namespace ROL2 {

template<class T>
constexpr bool is_enum_v = std::is_enum<T>::value;

template<class ENumType, class ReturnType>
using enable_if_enum_t = std::enable_if_t<is_enum_v<ENumType>,ReturnType>;

template<class ENumType>
constexpr auto enum_size_v = static_cast<std::underlying_type_t<ENumType>>(ENumType::Last);

template<typename T>
using member_Type_t = typename T::Type;

template<typename T>
using member_Flag_t = typename T::Flag;


namespace detail {

// Detect if T has a member type called Type that is an enum
template<typename T, typename=void>
struct has_enum_Type : std::false_type {};

template<typename T>
struct has_enum_Type<T,void_t<member_Type_t<T>> {
  static constexpr auto value = is_enum_v<member_Type_t<T>>;
};

// Detect if T has a member type called Flag that is an enum
template<typename T, typename=void>
struct has_enum_Flag : std::false_type {};

template<typename T>
struct has_enum_Flag<T,void_t<member_Flag_t<T>>> {
  static constexpr auto value = is_enum_v<member_Flag_t<T>>;
};

} // namespace detail

template<typename T>
constexpr bool has_enum_Type_v = detail::has_enum_Type<T>::value;


template<typename T>
constexpr bool has_enum_Flag_v = detail::has_enum_Type<T>::value;

template<typename object_type, typename return_type>
using enable_if_has_enum_Type_t = std::enable_if_t<has_enum_Type_v<object_type>,return_type>;

template<typename object_type, typename return_type>
using enable_if_has_enum_Flag_t = std::enable_if_t<has_enum_Type_v<object_type>,return_type>;


/** \brief Universal valid value function for all ROL2 enums
           All enums in ROL2 follow the same design

           enum class ENumType : SomeIntegralType {
             FirstValue = 0,
             ...
             Last 
           };
    
          The ENumType values always start at 0, are consecutive,
          and the largest value is denoted as Last. Values are 
          considered valid if they are at least 0 and less than Last.
          This function will throw an exception if the value is less than
          0 or greater than ENumType::Last
*/

template<class ENumType>
enable_if_enum_t<ENumType,bool>
is_valid_enum_value( const ENumType& e ) {
  using integral_type = std::underlying_type_t<ENumType>;
  constexpr auto first   = static_cast<integral_type>(0);
  constexpr auto last    = static_cast<integral_type>(ENumType::Last);
            auto current = static_cast<integral_type>(e);
  ROL_TEST_FOR_EXCEPTION( (current < first) || (current > last) ),
                          std::out_of_range, "The enum value is out of range" );
  return (first <= current) && (current < last);
};

/** \brief Universal post-increment operator of ROL2 class enums
*/
template<typename ENumType>
enable_if_enum_t<ENumType,ENumType&&>
operator++ ( ENumType& e ) {
  using int_type = std::underlying_type_t<ENumType>;
  auto e_int = static_cast<int_type>(e);
  constexpr auto last = static_cast<int_type>(EType::Last);
  ROL_TEST_FOR_EXCEPTION( e_int >= last, std::out_of_range, 
                          "The enum value is out of range" );
  e = static_cast<ENumType>(++e_int);
  return e; 
}

/** \brief Universal pre-increment operator of ROL2 class enums
*/
template<typename ENumType>
enable_if_enum_t<ENumType,ENumType&&>
operator++ ( ENumType& e, int ) {
  using int_type = std::underlying_type_t<ENumType>;
  auto e_int = static_cast<int_type>(e);
  constexpr auto last = static_cast<int_type>(EType::Last);
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


#endif //ROL2_ENUMERATION_HPP

