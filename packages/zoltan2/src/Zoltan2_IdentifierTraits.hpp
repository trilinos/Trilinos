// @HEADER
// ***********************************************************************
//            copyright
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERTRAITS
#define _ZOLTAN2_IDENTIFIERTRAITS

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <Zoltan2_config.h>

/*! \file Zoltan2_IdentifierTraits.hpp
  \brief Defines basic traits for application supplied local and global IDs.

*/

namespace Z2
{

/*! \struct Z2::IdentifierTraits
  \brief This structure defines the basic traits needed for an identifier.

    If an identifier type is not instantiated here, the caller will get a
    compile time error message indicating that we need to add a new type.
*/


// Is the type a Teuchos global ordinal type

template <typename T>
bool isGlobalOrdinalType<T> {
  return Teuchos::OrdinalTraits<T>::hasMachineParameters;
}

// Make a compact string from an integer

template <typename T>
std::string hexString(const T &ival){ 
  ostringstream os;
  os << std::hex << ival;
  return os.str();
}

// More compact.  Will this work even though characters are unprintable?

template <typename T>
std::string byteString(const T &ival)
{
  std::string s;

  if (ival != 0){             // scheme that works if non-zero
    s.reserve(sizeof(T));
    for (int i=sizeof(T)-1; i>=0; i--){
      int movebits = i*8;
      char c = static_cast<char>((ival & (0xff << movebits)) >> movebits);
      if ((s.size() > 0) || (static_cast<int>(c)  != 0)){
        s.push_back(c);
      }
    }
  }
  else{
    s.reserve(sizeof(T)+1);  // unique value for 0
    for (int i=0; i <= sizeof(T); i++){
      s.push_back('-');
    }
  }

  return s;
}

template<typename T>
struct UndefinedIdentifierTraits
{
  static inline T notDefined() { return T::Identifer_Type_Needs_To_Be_Added(); }
};
	
template<typename T>
struct IdentifierTraits {

  // An integer for each id, for hashing to a value between 1 and n 
  static inline int hashCode(const T& ) { 
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  // A unique key for each id, for hashtable use
  static inline std::string key(const T&){
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }
		
  // The name of this Identifier type
  static inline std::string name() { 
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template<>
struct IdentifierTraits<char> {
  static inline int hashCode(const char &c) { return static_cast<int>(c);}
  static inline std::string key(const char&c){ return std:string(1, c); }
  static inline std::string name()     { return("char");}
};

template<>
struct IdentifierTraits<short int> {
  static inline int hashCode(const short int &a) { return static_cast<int>(a);}
  static inline std::string key(const short int &a){ return byteString(a); }
  static inline std::string name()          { return("short int");}
};

template<>
struct IdentifierTraits<int> {
  static inline int hashCode(const int &a) { return a; }
  static inline std::string key(const int &a){ return byteString(a); }
  static inline std::string name()    { return("int");}
};

template<>
struct IdentifierTraits<unsigned int> {
  static inline int hashCode(const unsigned int &a) { 
    return static_cast<int>(a);
  }
  static inline std::string key(const unsigned int &a){ return byteString(a); }
  static inline std::string name()             { return("unsigned int");}
};

template<>
struct IdentifierTraits<long int> {
  static inline int hashCode(const long int &a) { 
    unsigned int total=0;
    for (int i=0, bits=0; i < sizeof(long int); i+=2, bits += 8){
      total += static_cast<unsigned int>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline std::string key(const long int &a){ return byteString(a); }
  static inline std::string name()    { return("long int");}
};

template<>
struct IdentifierTraits<unsigned long int> {
  static inline int hashCode(const unsigned long int &a) { 
    return IdentiferTraits<long int>::hashCode(static_cast<const long int>(a));
  }
  static inline std::string key(const unsigned long int &a){ 
    return byteString(a); 
  }
  static inline std::string name()   { return("long unsigned int");}
};

template<>
struct IdentifierTraits<long long int> {
  static inline int hashCode(const long long int &a) { 
    unsigned int total=0;
    for (int i=0, bits=0; i < sizeof(long long int); i+=2, bits += 8){
      total += static_cast<unsigned int>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline std::string key(const long long int &a){ return byteString(a); }
  static inline std::string name()    { return("long long int");}
};

template<>
struct IdentifierTraits<unsigned long long int> {
  static inline int hashCode(const unsigned long long int &a) { 
    return IdentiferTraits<long long int>::hashCode(
      static_cast<const long long int>(a));
  }
  static inline std::string key(const unsigned long long int &a){ 
    return byteString(a); 
  }
  static inline std::string name()    { return("long long int");}
};

template<>
struct IdentifierTraits<void *> {
  static inline int hashCode(void * const p) {
    return IdentiferTraits<intptr_t>::hashCode(static_cast<const intptr_t>(p));
  }
  static inline std::string key(void * const p){ 
    return byteString(static_cast<const intptr_t>(p));
  }
  static inline std::string name()   { return "pointer"; } 
};

template<typename T1, typename T2>
struct IdentifierTraits<std::pair<T1, T2> > {
  static inline int hashCode(const std::pair<T1, T2> &p)  {
    return IdentiferTraits<T1>::hashCode(p.first) + 
      IdentiferTraits<T2>::hashCode(p.second);
  }

  static inline std::string key(const std::pair<T1, T2> &p)  {
    return byteString(p.first) + std::string(":") + hexString(p.second);
  }
    
  static inline std::string name() {
     return std::string("std::pair<") + 
       IdentifierTraits<T1>::name() + std::string(", ");
       IdentifierTraits<T2>::name() + std::string(">");
  }
};

#ifdef SERIALIZATION_SUPPORTS_VECTORS
template<typename T>
struct IdentifierTraits<std::vector<T > > {
  static inline int hashCode(const std::vector<T> &v){
    unsigned long total=0;
    for (size_t i=0; i < v.size(); i++){
      total += IdentifierTraits<T>::hashCode(v[i]); 
    }
    return static_cast<int>(total);
  }  
  static inline std::string key(const std::vector<T> &v){
    std::string s;
    if (v.size()){
      s = s + byteString(v[0]);
      for (size_t i=1; i < v.size(); i++){
        s = s + std::string(":") + byteString(v[i]);
      }
    }
    return std::string;
  }
  static inline std::string name(){
     return std::string("std::vector<") + 
       IdentifierTraits<T>::name() + std::string(">");
  }
};
#endif



#endif // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
