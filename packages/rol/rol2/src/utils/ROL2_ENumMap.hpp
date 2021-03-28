#pragma once
#ifndef ROL2_ENUMMAP_HPP
#define ROL2_ENUMMAP_HPP

// #include <typeinfo>

namespace ROL2 {

/** \class ROL2::EnumMap
 *  \brief Defines a standardized bijective dictionary type object for
 *         converting between ROL's option strings and their corresponding
           enum values
 */

template<typename EType>
class EnumMap {
public:
  EnumMap( std::initializer_list<std::string> ilist ) {
    using size_type = typename std::initializer_list<std::string>::size_type;
//    std::cout << "EnumType = " << typeid(EType).name() << std::endl;
//    std::cout << "EnumMap : size(ilist) == " << ilist.size() << ", Last = " << EType::Last << std::endl;
    assert( ilist.size() == static_cast<size_type>(EType::Last) );
    auto e = static_cast<EType>(0);

    for( const auto& s : ilist ) {
      s2e[remove_format(s)] = e;
      e2s[e] = s;
      e++;
    }
    e2s[EType::Last] = "Last Type";
    s2e["lasttype"]  = EType::Last;
  }

  inline const std::string& operator[] ( EType e ) const { 
    if(!is_valid(e)) {
      std::stringstream msg;
      msg << "Enum Key " << e << " is invalid" << std::endl;
      throw std::out_of_range(msg.str());
    }
    return e2s.at(e); 
  }

  inline EType operator[] ( std::string s ) const { 
    if(!is_valid(s)) {
      std::stringstream msg;
      msg << "String Key " << s << " is invalid" << std::endl;
      msg << "Valid Keys are" << std::endl;
      for( int i=0; i<=EType::Last; ++i ) 
      msg << "  " << (*this)[static_cast<EType>(i)] << " : " << i << std::endl;
      throw std::out_of_range(msg.str());
    }
    return s2e.at(s); 
  }

  inline bool is_valid( EType e ) const { return e2s.count(e) == 1; }

  inline bool is_valid( std::string s ) const { return s2e.count(remove_format(s)) == 1; }

private:

  static std::string remove_format( std::string s ) {
    s.erase(std::remove_if(s.begin(), s.end(),
      []( auto const& c ) -> bool { return !std::isalnum(c); } ), s.end());
    for( auto& c : s ) c = tolower(c);
    return s;
  }

  std::map<EType,std::string> e2s;
  std::map<std::string,EType> s2e;
};

//template<typename object_type>
//inline
//enable_if_has_enum_Type_t<object_type,member_Type_t<object_type>>
//stringToEnum( std::string s, const object_type& ) {
//  return object_type::type_dict[s];
//}
//
//template<typename object_type>
//inline
//enable_if_has_enum_Flag_t<object_type,member_Flag_t<object_type>>
//stringToEnum( std::string s, const object_type& ) {
//  return object_type::flag_dict[s];
//}

} // namespace ROL2
#endif // ROL2_ENUMMAP_HPP

