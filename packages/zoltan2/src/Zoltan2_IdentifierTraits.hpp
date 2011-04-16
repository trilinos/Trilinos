// @HEADER
// ***********************************************************************
//            copyright
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERTRAITS
#define _ZOLTAN2_IDENTIFIERTRAITS

#include <utility>
#include <iostream>
#include <Zoltan2_config.h>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_HashUtils.hpp>

/*! \file Zoltan2_IdentifierTraits.hpp
  \brief Defines basic traits for application supplied local and global IDs.

    Note: Do not try to write an Identifier value to an ostream.  
    Some Identifier data types can not be written to an ostream.
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
      
  static int hashCode(const T id) { 
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Compute a key which will be unique for each id.
      \param id an application's local or global identifier
      \result the key
   */
  
  static double key(const T id){
   return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Convert a key back to the identifier that generated it.
      \param the key
      \result the identifier that would have generated this key
   */

  static T keyToGid(const double x){
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

  static inline bool equal(const T a, const T b) {
    return UndefinedIdentifierTraits<T>::notDefined(); 
  }

  /*! \brief Determine if data type can be an application global ID.
      \result true if they can.
   */

  static inline bool is_valid_id_type() { return false; }
};

/*! \cond IndentifierTraits_definitions
 */

template<>
struct IdentifierTraits<char> {
  static inline int hashCode(const char c) { return static_cast<int>(c);}
  static inline double key(const char c){ return static_cast<double>(c); }
  static char keyToGid(const double key){ return static_cast<char>(key); }
  static inline std::string name()     { return("char");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket< char >(); }
  static inline bool equal(const char  a, const char  b) { return (a==b); }
  static inline bool is_valid_id_type() { return true; }
};

template<>
struct IdentifierTraits<short int> {
  static inline int hashCode(const short int  a) { return static_cast<int>(a);}
  static inline double key(const short int a){ 
    return static_cast<double>(a); }
  static short int keyToGid(const double key){ 
    return static_cast<short int>(key);}
  static inline std::string name()          { return("short int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<short int>(); }
  static inline bool equal(const short int a, const short int b) { 
    return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

template<>
struct IdentifierTraits<int> {
  static inline int hashCode(const int a) { return a; }
  static inline double key(const int a){ return static_cast<double>(a); }
  static int keyToGid(const double key){ return static_cast<int>(key); }
  static inline std::string name()    { return("int");}
  static inline bool isHashKeyType() { return true; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<int>(); }
  static inline bool equal(const  int a, const  int b) { 
    return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

template<>
struct IdentifierTraits<unsigned int> {
  static inline int hashCode(const unsigned int a) { 
    return static_cast<int>(a); }
  static inline double key(const unsigned int a){ 
    return static_cast<double>(a); }
  static unsigned int keyToGid(const double key){  
    return static_cast<unsigned>(key);}
  static inline std::string name()             { return("unsigned int");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<unsigned int>(); }
  static inline bool equal(const unsigned int a, const unsigned int b) { 
    return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

template<>
struct IdentifierTraits<long> {
  static inline int hashCode(const long a) { 
    unsigned int total=0;
    for (unsigned int i=0, bits=0; i < sizeof(long); i++, bits += 8){
      total += static_cast<unsigned int>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const long a){ return static_cast<double>(a); }
  static long keyToGid(const double key){ return static_cast<long>(key); }
  static inline std::string name()    { return("long");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<long>(); }
  static inline bool equal(const long a, const long b) { 
    return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

template<>
struct IdentifierTraits<unsigned long> {
  static inline int hashCode(const unsigned long a) { 
    return IdentifierTraits<long>::hashCode(static_cast<unsigned long>(a)); }
  static inline double key(const unsigned long a){ 
    return static_cast<double>(a);}
  static unsigned long keyToGid(const double key){ 
    return static_cast<unsigned long>(key);}
  static inline std::string name()   { return("unsigned long");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<unsigned long>(); }
  static inline bool equal( const unsigned long a, 
    const unsigned long b) { return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

// TODO ifdef LONG LONG

template<>
struct IdentifierTraits<long long> {
  static inline int hashCode(const long long a) { 
    unsigned int total=0;
    for (unsigned int i=0, bits=0; i < sizeof(long long); i++, bits += 8){
      total += static_cast<unsigned int>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const long long a){ 
    return static_cast<double>(a);}
  static long long keyToGid(const double key){ 
    return static_cast<long long>(key);}
  static inline std::string name()    { return("long long");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return true; }
  static inline bool isPacketType() { return isPacket<long long>(); }
  static inline bool equal( const long long a, const long long b) { 
    return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

template<>
struct IdentifierTraits<unsigned long long > {
  static inline int hashCode(const unsigned long long a) { 
    return IdentifierTraits<long long>::hashCode(static_cast<long long>(a));}
  static inline double key(const unsigned long long a){ 
    return static_cast<double>(a); }
  static unsigned long long keyToGid(const double key){ 
    return static_cast<unsigned long long>(key);}
  static inline std::string name()    { return("unsigned long long");}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinalType() { return false; }
  static inline bool isPacketType() { 
    return isPacket<unsigned long long>(); }
  static inline bool equal( const unsigned long long a, 
    const unsigned long long b) { return (a==b) ; }
  static inline bool is_valid_id_type() { return true; }
};

template<typename T1, typename T2>
struct IdentifierTraits<std::pair<T1, T2> > {
  static inline int hashCode(const std::pair<T1, T2> p)  {
    return IdentifierTraits<T1>::hashCode(p.first) + 
      IdentifierTraits<T2>::hashCode(p.second);
  }

  static inline double key(const std::pair<T1, T2> p)  {
    int nbits = sizeof(T1)*8;
    intmax_t part1(p.first);
    intmax_t part2(p.second);
    return static_cast<double>((part2 << nbits) | part1);
  }

  static std::pair<T1, T2> keyToGid(const double key){
    intmax_t t1_size_mask = 0xff;
    for (unsigned int n=1; n < sizeof(T1); n++){
      t1_size_mask |= (0xff << n);
    }
    intmax_t ikey = static_cast<intmax_t>(key);
    T1 part1 = static_cast<T1>(ikey & t1_size_mask);
    T2 part2 = static_cast<T2>((ikey & ~t1_size_mask) >> sizeof(T1)*8);
    return std::pair<T1,T2>(part1, part2);
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
  static inline bool equal( const std::pair<T1, T2> a, 
    const std::pair<T1, T2> b) { 
    return ((a.first==b.first) && (a.second==b.second)); }
  static inline bool is_valid_id_type() { 
    return (sizeof(T1)+sizeof(T2) <= sizeof(intmax_t)); }
};

/*! \endcond
 */

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
