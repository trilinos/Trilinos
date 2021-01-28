#pragma once
#ifndef ROL2_PARAMETERLIST_HPP
#define ROL2_PARAMETERLIST_HPP

#include <map>
#include <string>
#include <memory>

namespace ROL2 {
class ParameterList {
public:
  using std::string;
  template<typename value_type> using dict = std::map<string,value_type>;

  //-----------------------------------------------------------------------
  // Fundamental Type Access

  inline void set( string key, bool   value ) { bools_[key]   = value; }
  inline void set( string key, int    value ) { ints_[key]    = value; }
  inline void set( string key, double value ) { doubles_[key] = value; }
  inline void set( string key, string value ) { strings_[key] = value; }

  // These throw std::out_of_range exceptions if key not in dict
  inline bool   get_bool(   string key ) const { return bools_.at(key);   }
  inline int    get_int(    string key ) const { return ints_.at(key);    }
  inline double get_double( string key ) const { return doubles_.at(key); }
  inline string get_string( string key ) const { return strings_.at(key); }

  bool get( string key, bool default_value ) { 
    if( !bools_.count(key) ) set( key, default_value );
    return get_bool(key);
  }

  int get( string key, int default_value ) { 
    if( !ints_.count(key) ) set( key, default_value );
    return get_ints(key);
  }

  double get( string key, double default_value ) { 
    if( !doubles_.count(key) ) set( key, default_value );
    return get_double(key);
  }

  string get( string key, string default_value ) { 
    if( !strings_.count(key) ) set( key, default_value );
    return get_string(key);
  }

  //-----------------------------------------------------------------------
  // Sublist Access

  ParameterList& sublist( string key, bool mustAlreadyExist=false ) {
    if( mustAlreadyExist ) return *(sublists_.at(key));
    else {   
      if( !sublist_.count(key) ) sublists_[key] = std::make_unique<ParameterList>();
      return *(sublists_[key]);
    }
  }

  const ParameterList& sublist( string key ) const {
    return *(sublists_.at(key));
  }

private:

  dict<bool>                           bools_;
  dict<int>                            ints_;
  dict<double>                         doubles_;
  dict<string>                         strings_;
  dict<std::unique_ptr<ParameterList>> sublists_;

}; // ParameterList
} // namespace ROL

#endif // ROL2_PARAMETERLIST_HPP


