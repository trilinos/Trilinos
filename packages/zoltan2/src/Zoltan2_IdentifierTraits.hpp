// @HEADER
// ***********************************************************************
//            copyright
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERTRAITS
#define _ZOLTAN2_IDENTIFIERTRAITS

#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <Zoltan2_config.h>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_HashUtils.hpp>

/*! \file Zoltan2_IdentifierTraits.hpp
  \brief Defines basic traits for application supplied local and global IDs.

    Note: Do not try to write an Identifier value to an ostream.  
    Instead write the string IdentifierMap<>::key().  Some Identifier
    data types can not be written to an ostream.
*/

/*! \namespace Z2
  \brief Internal Zoltan2 namespace.

  This namespace contains the objects and functions that are not
  part of the Zoltan2 user interface.
*/

namespace Z2
{

/*! 
    \brief Determine whether id type is a Teuchos communication Packet type.
    \tparam T data type for identifier
    \result true if it is a Teuchos Packet type
 */

template <typename T>
bool isPacket() {
  return Teuchos::SerializationTraits<int, T>::supportsDirectSerialization;
}


 /* \brief Represent a numeric value as a hex string.

    The goal was to come up with a way to uniquely represent identifiers
    compactly in a string.  Some values would require fewer bytes when 
    represented as a hex string.

    \tparam T the numeric data type to convert to the hex string
    \result the hex string (e.g. 17 -> "11")
 */

template <typename T>
std::string keyToHexString(const T &ival){ 
  std::ostringstream os;
  os << std::hex << ival;
  return os.str();
}

template <typename T>
T hexStringToKey(const std::string s){
  T val;
  std::istringstream is(s);
  is >> std::hex >> val;
  return val;
}

/*! 
    \brief Represent a numeric value as a byte string.

    The goal was to come up with a way to represent identifiers
    compactly in a string.  This scheme is more compact than
    hexString, but includes unprintable characters in the string.
 
    This is a problem when debugging, and Hashtable error messages
    try to print out keys.

    \tparam T the numeric data type to convert to the string
    \param ival the number to convert to a string
    \result the string representation.
 */

template <typename T>
std::string keyToByteString(const T &ival)
{
  std::string s;

  if (ival != 0){             // scheme that works if non-zero
    s.reserve(sizeof(T));
    for (unsigned int i=sizeof(T); i>0; i--){
      int movebits = (i-1)*8;
      char c = static_cast<char>((ival & (0xff << movebits)) >> movebits);
      if ((s.size() > 0) || (static_cast<int>(c)  != 0)){
        s.push_back(c);
      }
    }
  }
  else{
    s = std::string("zero");  // unique value for 0
    for (unsigned int i=s.size(); i <= sizeof(T); i++){
      s.push_back('!');
    }
  }

  return s;
}

template <typename T>
T byteStringToKey(std::string s)
{
  // special case of zero

  if (s.size() == (sizeof(T) + 1)){
    return 0;
  }

  // everything else
  
  T val = 0;
  int shift = 0;
  std::string::iterator next = s.end();
  std::string::iterator start = s.begin();

  while (next != start){
    int ic = static_cast<int>(*next);
    val = val | (ic << shift);
    shift += 8;
    next --;
  }
  return val;
}

template <typename T>
std::string keyToString(const T &ival){ 
#ifdef Z2_PRINTABLE_HASH_TABLE_STRING
  return keyToHexString<T>(ival);
#else
  return keyToByteString<T>(ival);
#endif
}
template <typename T>
T stringToKey(const std::string s){
#ifdef Z2_PRINTABLE_HASH_TABLE_STRING
  return hexStringToKey<T>(s);
#else
  return byteStringToKey<T>(s);
#endif
}

/*! \struct UndefinedIdentifierTraits
    \brief Resolves at compile time the use of invalid id types.
*/

template<typename T>
struct UndefinedIdentifierTraits
{ 
static inline T notDefined() { return T::Identifier_Type_Needs_To_Be_Added(); } };


/*! \struct IdentifierTraits
    \brief The functions to be defined for each valid id type.
    \tparam  T the local or global identifier data type
*/

template<typename T>
struct IdentifierTraits {

  /*! \brief Compute an integer hash code for the id.
      \param id an application's local or global identifier
      \result the integer code, need not be unique for each id
   */
      
  static int hashCode(const T& id) { 
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Compute a key which will be unique for each id.
      \param id an application's local or global identifier
      \result the key, as a std::string
   */
  
  static std::string key(const T& id){
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Convert a key back to the identifier that generated it.
      \param key a string created with a call to IdentifierTraits::key()
      \result the identifier that would have generated this key
   */

  static T keyToGid(const std::string key){
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }
		
  /*! \brief The name of the identifier data type.
      \result The name as a std::string.
   */

  static inline std::string name() { 
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Inform caller if the data type can be used by Teuchos::hashCode.
      \result true if the data type can be used with Teuchos::hashCode.
   */

  static inline bool isHashKeyType() {
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Inform caller if the data type is a Teuchos GlobalOrdinal
      \result true if it is a Teuchos GlobalOrdinal
   */

  static inline bool isGlobalOrdinalType() {
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Inform caller if the data type is a Teuchos Packet type 
      \result true if it is a Teuchos Packet type
   */

  static inline bool isPacketType() {
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Determine if two identifiers are the same.
      \result true if they are the same.
   */

  static inline bool equal(const T &a, const T &b) {
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }
};

/*! \cond IndentifierTraits_definitions
 */

template<>
struct IdentifierTraits<char> {
  static inline int hashCode(const char &c) { return static_cast<int>(c);}
  static inline std::string key(const char&c){ return keyToString<char>(c); }
  static char keyToGid(const std::string s){ return stringToKey<char>(s); }
  static inline std::string name()     { return("char");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket< char >(); }
  static inline bool equal(const char &a, const char &b) { return (a==b); }
};

template<>
struct IdentifierTraits<short int> {
  static inline int hashCode(const short int &a) { return static_cast<int>(a);}
  static inline std::string key(const short int &a){ return keyToString(a); }
  static short int keyToGid(const std::string s){ 
    return stringToKey<short int>(s);}
  static inline std::string name()          { return("short int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<short int>(); }
  static inline bool equal(const short int &a, const short int &b) { 
    return (a==b) ; }
};

template<>
struct IdentifierTraits<int> {
  static inline int hashCode(const int &a) { return a; }
  static inline std::string key(const int &a){ return keyToString(a); }
  static int keyToGid(const std::string s){ 
    return stringToKey<int>(s);}
  static inline std::string name()    { return("int");}
  static inline bool isHashKeyType() { return true; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<int>(); }
  static inline bool equal(const  int &a, const  int &b) { 
    return (a==b) ; }
};

template<>
struct IdentifierTraits<unsigned int> {
  static inline int hashCode(const unsigned int &a) { 
    return static_cast<int>(a); }
  static inline std::string key(const unsigned int &a){ return keyToString(a); }
  static unsigned int keyToGid(const std::string s){ 
    return stringToKey<unsigned int>(s);}
  static inline std::string name()             { return("unsigned int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<unsigned int>(); }
  static inline bool equal(const unsigned int &a, const unsigned int &b) { 
    return (a==b) ; }
};

template<>
struct IdentifierTraits<long int> {
  static inline int hashCode(const long int &a) { 
    unsigned int total=0;
    for (unsigned int i=0, bits=0; i < sizeof(long int); i+=2, bits += 8){
      total += static_cast<unsigned int>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline std::string key(const long int &a){ return keyToString(a); }
  static long int keyToGid(const std::string s){ 
    return stringToKey<long int>(s);}
  static inline std::string name()    { return("long int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<long int>(); }
  static inline bool equal(const long int &a, const long int &b) { 
    return (a==b) ; }
};

template<>
struct IdentifierTraits<unsigned long int> {
  static inline int hashCode(const unsigned long int &a) { 
    return 
      IdentifierTraits<long int>::hashCode(static_cast<const long int>(a)); }
  static inline std::string key(const unsigned long int &a){ 
    return keyToString(a); }
  static unsigned long int keyToGid(const std::string s){ 
    return stringToKey<unsigned long int>(s);}
  static inline std::string name()   { return("long unsigned int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<unsigned long int>(); }
  static inline bool equal( const unsigned long int &a, 
    const unsigned long int &b) { return (a==b) ; }
};

template<>
struct IdentifierTraits<long long int> {
  static inline int hashCode(const long long int &a) { 
    unsigned int total=0;
    for (unsigned int i=0, bits=0; i < sizeof(long long int); i+=2, bits += 8){
      total += static_cast<unsigned int>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline std::string key(const long long int &a){ return keyToString(a); }
  static long long int keyToGid(const std::string s){ 
    return stringToKey<long long int>(s);}
  static inline std::string name()    { return("long long int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<long long int>(); }
  static inline bool equal( const long long int &a, const long long int &b) { 
    return (a==b) ; }
};

template<>
struct IdentifierTraits<unsigned long long int> {
  static inline int hashCode(const unsigned long long int &a) { 
    return IdentifierTraits<long long int>::hashCode(
      static_cast<const long long int>(a)); }
  static inline std::string key(const unsigned long long int &a){ 
    return keyToString(a); }
  static unsigned long long int keyToGid(const std::string s){ 
    return stringToKey<unsigned long long int>(s);}
  static inline std::string name()    { return("long long int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return false; }
  static inline bool isPacketType() { 
    return isPacket<unsigned long long int>(); }
  static inline bool equal( const unsigned long long int &a, 
    const unsigned long long int &b) { return (a==b) ; }
};

template<typename T1, typename T2>
struct IdentifierTraits<std::pair<T1, T2> > {
  static inline int hashCode(const std::pair<T1, T2> &p)  {
    return IdentifierTraits<T1>::hashCode(p.first) + 
      IdentifierTraits<T2>::hashCode(p.second);
  }

  static inline std::string key(const std::pair<T1, T2> &p)  {
    return keyToString(p.first) + std::string(":") + keyToString(p.second);
  }

  static std::pair<T1, T2> keyToGid(const std::string s){
    std::string::size_type pos = s.find(":");
    T1 first = stringToKey<T1>(s.substr(0,pos));
    T2 second = stringToKey<T2>(s.substr(pos+1));
    std::pair<T1, T2> p(first, second);
    return p;
  }
    
  static inline std::string name() {
     std::string s("std::pair<");
     s += IdentifierTraits<T1>::name();
     s += std::string(", ");
     s += IdentifierTraits<T2>::name();
     s += std::string(">");
     return s;
  }
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return false; }
  static inline bool isPacketType() { return isPacket<std::pair<T1, T2> >(); }
  static inline bool equal( const std::pair<T1, T2> &a, 
    const std::pair<T1, T2> &b) { 
    return ((a.first==b.first) && (a.second==b.second)); }
};

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
      s = s + keyToString(v[0]);
      for (size_t i=1; i < v.size(); i++){
        s = s + std::string(":") + keyToString(v[i]);
      }
    }
    return s;
  }

  static std::vector<T> keyToGid(const std::string s){
    std::string scopy(s);
    int num=0;
    std::string::size_type k = scopy.find(":");
    while (k != std::string::npos){
      num++;
      scopy = scopy.substr(k+1);
      k = scopy.find(":");
    }
    std::vector<T> v;
    v.reserve(num+1);
    scopy = s;
    k = scopy.find(":");
    while (k != std::string::npos){
      v.push_back(stringToKey<T>(scopy.substr(0, k)));
      scopy = scopy.substr(k+1);
      k = scopy.find(":");
    }
    v.push_back(stringToKey<T>(scopy.substr(0, k)));
    return v;
  }
  static inline std::string name(){
     return std::string("std::vector<") + 
       IdentifierTraits<T>::name() + std::string(">");
  }
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return false; }
  static inline bool isPacketType() { return isPacket<std::vector<T> >(); }

  static inline bool equal(const std::vector<T> &a, const std::vector<T> &b){
    if (a.size() != b.size()) return false;
    for (unsigned int i=0; i < a.size(); i++)
      if (a[i] != b[i]) return false;
    return true;
  }
};

/*! \endcond
 */

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
