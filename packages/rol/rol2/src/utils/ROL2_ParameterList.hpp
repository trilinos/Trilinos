#pragma once
#ifndef ROL2_PARAMETERLIST_HPP
#define ROL2_PARAMETERLIST_HPP

namespace ROL2 {

namespace detail {

template<typename value_type>
std::ostream& operator << ( std::ostream& os, std::vector<value_type>& v ) {
  os << "[ ";
  for( std::size_t k=0; k<v.size()-1; ++k ) os << v[k] << ", ";
  os << v.back() << " ]" << std::endl;
  return os;
}


template<typename value_type>
class PList {
public:

  inline void set( std::string key, value_type value ) {
    values_[key] = value;
  }
 
  inline value_type get( std::string key, value_type value ) {
    if( values_.count(key) ) values_[key] = value;
    return values_[key];
  }

  virtual int get_level() = 0;// { return 0; }

  friend std::ostream& operator << ( std::ostream& os, PList& plist ) {
    for( auto it=plist.values_.begin(); it != plist.values_.end(); ++it )
      os << std::string(2*(plist.get_level()),' ') << 
            it->first << " : " << it->second << std::endl;
    return os;
  }

private:
  std::map<std::string,value_type> values_;  
}; // class detail::PList



template<typename...ValueTypes>
class ParameterList : public PList<ValueTypes>... {
public:

  ParameterList( int level = 0 ) : level_(level) {}

  template<typename value_type>
  inline void set( std::string key, value_type value ) {
    static_assert( disjunction_v<is_convertible<value_type,ValueTypes>...>,"" );
    static_cast<PList<value_type>&>(*this).set( key, value );
  }

  inline void set( std::string key, const char* value ) {
    set(key,std::string(value));
  }

  template<typename value_type>
  inline value_type get( std::string key, value_type value ) {
    static_assert( disjunction_v<is_convertible<value_type,ValueTypes>...>,"" );
    return PList<value_type>::get( key, value );
  }

  inline std::string get( std::string key, const char* value ) {
    return get(key,std::string(value));
  }

private:
  int level_ = 0;

}; // class detail::ParameterList

template<typename...ValueTypes>
void display( std::ostream& os, ParameterList<ValueTypes...>& parlist ) {
  using expander = int[];
  (void)expander{0, (void( os << static_cast<PList<ValueTypes>&>(parlist) ),0) ... };
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

  ParameterList( int level = 0 ) : level_(level) {}

  ParameterList& sublist( std::string key, bool mustAlreadyExist=false );

  int get_level() override { return level_; }

  friend std::ostream& operator << ( std::ostream&  os, 
                                     ParameterList& parlist ) {
    detail::display( os , parlist );
    os << std::endl;
    for( auto it=parlist.sublists_.begin(); it != parlist.sublists_.end(); ++it )
      os << std::string(2*parlist.get_level(),' ') <<
           it->first << ":\n" << *(it->second);
    return os;
  }

private:
  std::map<std::string,std::unique_ptr<ParameterList>> sublists_;
  int level_ = 0;  

}; // ParameterList

std::ostream& operator << ( std::ostream& os, ParameterList& parlist );

} // namespace ROL2

#endif // ROL2_PARAMETERLIST_HPP


