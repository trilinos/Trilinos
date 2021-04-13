#pragma once
#ifndef ROL_PARLIST_HPP
#define ROL_PARLIST_HPP

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
std::ostream& operator << ( std::ostream& os, std::vector<value_type>& v ) {
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
    if( !values_.count(key) ) values_[key] = value;
    return values_.at(key);
  }

  virtual int get_level() = 0;

  friend std::ostream& operator << ( std::ostream& os, 
                                     PList&        plist ) {
    for( auto it=plist.values_.begin(); it != plist.values_.end(); ++it )
      os << std::string(2*(plist.get_level()),' ') << 
            it->first << " : " << it->second << std::endl;
    return os;
  }

private:

  std::map<key_type,value_type> values_;  

}; // class detail::PList



template<typename...ValueTypes>
class ParList : public PList<ValueTypes>... {
public:

  using key_type = std::string;

  // True if type can be converted to any of the ParList value types
  template<typename value_type>
  static constexpr bool is_valid_type_v = disjunction_v<std::is_convertible<value_type,ValueTypes>...>;

  ParList( int level = 0 ) : level_(level) {}

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

  template<typename value_type>
  inline value_type get( key_type   key, 
                         value_type value ) {
    static_assert( is_valid_type_v<value_type>, "" );
    return PList<value_type>::get( key, value );
  }

  inline std::string get( key_type    key, 
                          const char* value ) {
    return get(key,key_type(value));
  }

private:
  int level_ = 0;

}; // class detail::ParList

template<typename...ValueTypes>
void display( std::ostream& os, 
              ParList<ValueTypes...>& parlist ) {
  using expander = int[];
  (void)expander{0, (void( os << static_cast<PList<ValueTypes>&>(parlist) ),0) ... };
}

} // namespace detail 


class ParList : public detail::ParList<bool, 
                                       int,
                                       double,
                                       std::string,
                                       std::vector<int>,
                                       std::vector<double>,
                                       std::vector<std::string>> {
public:

  using key_type = std::string;

  ParList( int level = 0 ) : level_(level) {}

  ParList& sublist( key_type key, 
                    bool     mustAlreadyExist=false );

  int get_level() override { return level_; }

  friend std::ostream& operator << ( std::ostream&  os, 
                                     ParList&       parlist ) {
    detail::display( os , parlist );
    os << std::endl;
    for( auto it=parlist.sublists_.begin(); it != parlist.sublists_.end(); ++it )
      os << std::string(2*parlist.get_level(),' ') <<
           it->first << ":\n" << *(it->second);
    return os;
  }

private:

  std::map<key_type,std::unique_ptr<ParList>> sublists_;
  int level_ = 0;  

}; // ParList

std::ostream& operator << ( std::ostream& os, ParList& parlist );

} // namespace ROL

#include "ROL_ParList.cpp"

#endif // ROL_PARLIST_HPP


