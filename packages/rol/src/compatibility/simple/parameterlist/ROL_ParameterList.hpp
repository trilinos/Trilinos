// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef ROL_PARAMETERLIST_HPP
#define ROL_PARAMETERLIST_HPP

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <type_traits>

#include "ROL_Stream.hpp"

namespace ROL {

#if __cplusplus >= 201703L 
template<class...B>
inline constexpr bool disjunction_v = std::disjunction<B...>::value;
#else
template<class...> struct disjunction : std::false_type { };
template<class B1> struct disjunction<B1> : B1 { };
template<class B1, class... Bn>
struct disjunction<B1, Bn...> 
    : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>>  {};
template<class...B>
constexpr bool disjunction_v = disjunction<B...>::value;
#endif

namespace detail {

template<typename value_type>
std::ostream& operator << ( std::ostream& os, const std::vector<value_type>& v ) {
  os << "[ ";
  for( std::size_t k=0; k<v.size()-1; ++k ) os << v[k] << ", ";
  os << v.back() << " ]" << std::endl;
  return os;
}


template<typename ValueType> 
class PList {
public:
 
  using key_type   = std::string;
  using value_type = ValueType;
  
  inline void set( key_type   key, 
                   value_type value ) {
    values_[key] = value;
  }
 
  inline value_type get( key_type   key, 
                         value_type value ) {
    if( values_.find(key) == values_.end() ) values_[key] = value;
    return values_.at(key);
  }

  inline value_type get( key_type key ) const {
    return values_.at(key);
  }

  virtual int get_level() const = 0;

  friend std::ostream& operator << ( std::ostream& os, 
                                     const PList& plist ) {
    for( auto it=plist.values_.begin(); 
              it != plist.values_.end(); ++it )
      os << std::string(2*(plist.get_level()),' ') << 
            it->first << " : " << it->second << std::endl;
    return os;
  }

private:
  const std::map<key_type,value_type>& get_values() const { return values_; }
  std::map<key_type,value_type> values_;  

}; // class detail::PList



template<typename...ValueTypes>
class ParameterList : public PList<ValueTypes>... {
public:

  using key_type = std::string;

  // True if type can be converted to any of the ParameterList value types
  template<typename value_type>
  static constexpr bool is_valid_type_v = disjunction_v<std::is_convertible<value_type,ValueTypes>...>;

  ParameterList( int level = 0 ) : level_(level) {}

  template<typename value_type>
  inline void set( key_type   key, 
                   value_type value ) {
    static_assert( is_valid_type_v<value_type>, "" );
    static_cast<PList<value_type>&>(*this).set( key, value );
  }

  inline void set( key_type    key, 
                   const char* value ) {
    set(key,key_type(value));
  }

  inline std::string get( key_type    key, 
                          const char* value ) {
    return get(key,std::string(value));
  }

  template<typename value_type>
  inline 
  std::enable_if_t<is_valid_type_v<value_type>,value_type> 
  get( key_type key ) const {
    //static_assert( is_valid_type_v<value_type>, "" );
    return PList<value_type>::get(key);
  }

  template<typename value_type>
  inline value_type get( key_type key, value_type value ) {
    //static_assert( is_valid_type_v<value_type>, "" );
    return PList<value_type>::get(key,value);
  }


private:
  int level_ = 0;

}; // class detail::ParameterList

template<typename...ValueTypes>
inline void display( std::ostream& os, 
                     const ParameterList<ValueTypes...>& parlist ) {
  using expander = int[];
  (void)expander{0, (void( os << static_cast<const PList<ValueTypes>&>(parlist) ),0) ... };
}

} // namespace detail 


class ParameterList : public detail::ParameterList<bool, 
                                                   int,
                                                   double,
                                                   std::string,
                                                   std::vector<int>,
                                                   std::vector<double>,
                                                   std::vector<std::string>> {
public:

  using key_type = std::string;

  ParameterList( int level = 0 ) : level_(level) {}

  ParameterList& sublist( key_type key, 
                          bool     mustAlreadyExist=false );

//  template<typename value_type>
//  inline value_type get( key_type key ) {
//    static_assert( is_valid_type_v<value_type>, "" );
//    return PList<value_type>::get( key );
//  }
//
//  template<typename value_type>
//  inline value_type get( key_type key, value_type default_value ) {
//    static_assert( is_valid_type_v<value_type>, "" );
//    return PList<value_type>::get( key,default_value );
//  }

  inline bool isSublist( key_type key ) {
    return sublists_.find(key) != sublists_.end();
  }

  int get_level() const override { return level_; }

  friend std::ostream& operator << ( std::ostream&  os, 
                                     const ParameterList& parlist ) {
    detail::display( os , parlist );
    os << std::endl;
    for( auto it=parlist.sublists_.begin(); it != parlist.sublists_.end(); ++it )
      os << std::string(2*parlist.get_level(),' ') <<
           it->first << ":\n" << *(it->second);
    return os;
  }

private:

  std::map<key_type,std::shared_ptr<ParameterList>> sublists_;
  int level_ = 0;  

}; // ParameterList

template<typename value_type>
inline std::vector<value_type>
getArrayFromStringParameter( const ParameterList& parlist, std::string key ) {
  return parlist.get<std::vector<value_type>>(key);
}


std::ostream& operator << ( std::ostream& os, const ParameterList& parlist );

} // namespace ROL

#include "ROL_ParameterList.cpp"

#endif // ROL_PARAMETERLIST_HPP


