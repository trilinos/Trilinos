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
#include <sstream>
#include <string>
#include <stdint.h>
#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_HashUtils.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_AlltoAll.hpp>

using Teuchos::SerializationTraits;

/*! \file Zoltan2_IdentifierTraits.hpp
  \brief Defines basic traits for application supplied local and global IDs.

  The data types permitted for global identifiers for Zoltan2 callers include
  some that are not represented in Teuchos::OrdinalTraits.  A common case is
  when a matrix nonzero is represented as an (i,j) pair.

  Zoltan2 uses the IdentifierTraits structures to manage the application
  supplied identifiers.
*/

/*! \namespace Zoltan2
  \brief Internal Zoltan2 namespace.

  This namespace contains the internal symbols that are not
  part of the Zoltan2 user interface.
*/
namespace Zoltan2
{


template <typename T>
  std::pair<T, T> z2LocalMinMax(const T *val, size_t n)
{
  T min = val[0], max = val[0];
  for (size_t i=1; i < n; i++){
    if (val[i] < min) min = val[i];
    else if (val[i] > max) max = val[i];
  }
  return std::pair<T,T>(min,max);
}

template <typename T>
  std::pair<T, T> z2GlobalMinMax(T minval, T maxval, const Comm<int> &comm)
{
  // We don't use Teuchos::reduceAll because T may not be a packet type.
  int nprocs = comm.getSize();
  Array<T> sendBuf(2*nprocs);
  for (int i=0; i < 2*nprocs; i+=2){
    sendBuf[i] = minval;
    sendBuf[i+1] = maxval;
  }
  ArrayRCP<T> recvBuf;
  Environment env;
  try{
    ArrayView<const T> sendView = sendBuf(); 
    AlltoAll<T, int>(comm, env, sendView, int(2), recvBuf);
  }
  catch(std::exception &e){
    Z2_THROW_ZOLTAN2_ERROR(env, e);
  }

  T *minPtr = recvBuf.getRawPtr();
  T *maxPtr = minPtr + 1;
  T min = *minPtr;
  T max = *maxPtr;
  
  for (int i=1; i < comm.getSize(); i++){
    minPtr += 2;
    maxPtr += 2;
    if (*minPtr < min) min = *minPtr;
    if (*maxPtr > max) max = *maxPtr;
  }

  return std::pair<T,T>(min, max);
}


template <typename T>
  bool z2AreConsecutive(const T *val, size_t n)
{
  if (val[n-1] - val[0] + 1 != n)
    return false;

  T nextval = val[0]+1;
  for (size_t i=1; i < n-1; i++){
    if (val[i] != nextval)
      return false;
    nextval++;
  }
  return true;
}

template <typename T>
  std::string stringifyOrdinal(T ordinal)
{
  std::ostringstream oss;
  oss << ordinal;
  return oss.str();
}

/*! 
    \brief Structure to catch invalid Indentifier types.
 */
template<typename T>
struct UndefIdTraits
{ 
static inline T invalid() { return T::Error(); } 
};


/*! \struct IdentifierTraits
    \brief The functions to be defined for each valid id type.
    \param  T the identifier data type
*/

template<typename T>
struct IdentifierTraits {

  /*! \brief Compute an integer hash code for the id.
      \param the id
      \result the integer code, need not be unique for each id
   */
  static int hashCode(const T id) { 
   return UndefIdTraits<int>::invalid(); 
  }

  /*! \brief Compute a key which will be unique for each id.
      \param the id
      \result the key
   */
  static double key(const T id){
   return UndefIdTraits<double>::invalid(); 
  }

  /*! \brief Convert a key back to the identifier that generated it.
      \param the key
      \result the identifier that would have generated this key
   */
  static T keyToGid(const double x){
   return UndefIdTraits<T>::invalid(); 
  }
		
  /*! \brief The name of the identifier data type.
      \result The name
   */
  static inline std::string name() { 
    return UndefIdTraits<std::string>::invalid(); 
  }
		
  /*! \brief A string representing the value.
      \result The  string
   */
  static std::string stringify(T val) { 
    return UndefIdTraits<std::string>::invalid(); 
  }

  /*! \brief Determine whether the data type can be used by Teuchos::hashCode.
      \result true if the data type can be used with Teuchos::hashCode
   */
  static inline bool isHashKeyType() {
    return UndefIdTraits<bool>::invalid(); 
  }

  /*! \brief Determine whether the data type can be a Teuchos Ordinal
      \result true if it can be a Teuchos Ordinal

      Only data types with definitions in Teuchos_OrdinalTraits.hpp
      can be used as Ordinals in Teuchos.
   */
  static inline bool isGlobalOrdinal() {
    return UndefIdTraits<bool>::invalid(); 
  }

  /*! \brief Cast to another type if possible
      \param x - an Identifier value
      \result y - x cast to another type if possible,
            otherwise create a run time error
   */
  static inline void castTo(const T x, char &y) {
      y = UndefIdTraits<char>::invalid();}

  static inline void castTo(const T x, short &y) {
      y = UndefIdTraits<short>::invalid();}

  static inline void castTo(const T x, int &y) {
      y = UndefIdTraits<int>::invalid();}

  static inline void castTo(const T x, unsigned &y) {
      y = UndefIdTraits<unsigned>::invalid();}

  static inline void castTo(const T x, long &y) {
      y = UndefIdTraits<long>::invalid();}

  static inline void castTo(const T x, unsigned long &y) {
      y = UndefIdTraits<unsigned long>::invalid();}

  static inline void castTo(const T x, long long &y) {
      y = UndefIdTraits<long long>::invalid();}

  static inline void castTo(const T x, unsigned long long &y) {
      y = UndefIdTraits<unsigned long long>::invalid();}

  template <typename T1, typename T2>
    static inline void castTo(const T x, std::pair<T1, T2> &y) { 
      y = UndefIdTraits<std::pair<T1, T2> >::invalid();}

  /*! 
      \brief Determine whether the data type can be used in 
                  Teuchos communication.
      \param T the data type 
      \result true if it can be a Teuchos Packet type

      Packet data types used in Teuchos_CommHelpers.hpp must have a 
     definition in SerializationTraits.
   */
  static inline bool isPacketType() {
    return UndefIdTraits<bool>::invalid(); 
  }

  /*! \brief Determine if two identifiers are the same type.
      \result true if they are the same.
   */
  static inline bool equal(const T a, const T b) {
    return UndefIdTraits<bool>::invalid(); 
  }

  /*! \brief Determine if a < b.
      \result true if a < b, false if !(a<b), and a
             compile error if T is not a data type that
             can be ordered
   */
  static inline bool lessThan(const T a, const T b) {
    return UndefIdTraits<bool>::invalid(); 
  }


  /*! \brief Determine if data type can be used by Zoltan2 caller an 
        application global ID.
      \result true if they can.
   */
  static inline bool is_valid_id_type() { return false; }

  /*! \brief Return the minimum and maximum of a list of values.
      \result A pair with the minimum value followed by the
        maximum value if T can be ordered, otherwise an error.
   */
  static std::pair<T, T> minMax(const T *values, size_t numValues) { 
    return UndefIdTraits<std::pair<T, T> >::invalid();
  }

  /*! \brief Return the global minimum and global maximum
      \param localMin - this process' local minimum
      \param localMax - this process' local maximum
      \param comm - A Teuchos::Comm for the reduction operation
      \result A pair with the global minimum value followed by the
        global maximum value if T can be ordered, otherwise an error.
   */
  static std::pair<T, T> globalMinMax(T localMin, T localMax, 
    const Comm<int> &comm) { 
      return UndefIdTraits<std::pair<T, T> >::invalid();
  }

  /*! \brief Determine if the values are increasing consecutive
      \param val - the list of values
      \param n - the number of values in the list
      \result Return true if the values are increasing consecutive,
            false if not, and if data type T can not be ordered,
            create a run time error
   */
  static bool areConsecutive(const T *val, size_t n){
    return UndefIdTraits<bool>::invalid();
  }

};

/*! \cond IndentifierTraits_definitions
 */

template<>
struct IdentifierTraits<char> {
  typedef unsigned long ulong;
  static inline int hashCode(const char c) { return static_cast<int>(c);}
  static inline double key(const char c){ return static_cast<double>(c); }
  static char keyToGid(const double key){ return static_cast<char>(key); }
  static inline std::string name()     { return("char");}
  static std::string stringify(char val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const char x, char &y) {y = x;}
  static inline void castTo(const char x, short &y) {y = short(x);}
  static inline void castTo(const char x, int &y) {y = int(x);}
  static inline void castTo(const char x, unsigned &y) {y = unsigned(x);}
  static inline void castTo(const char x, long &y) {y = long(x);}
  static inline void castTo(const char x, ulong &y) {y = ulong(x);}
#ifdef HAVE_LONG_LONG
  typedef unsigned long long ullong;
  typedef long long llong;
  static inline void castTo(const char x, llong &y) {y = llong(x);}
  static inline void castTo(const char x, ullong &y) {y = ullong(x);}
#endif
  template <typename T1, typename T2>
    static inline void castTo(const char x, std::pair<T1, T2> &y) { 
      throw std::runtime_error("invalid conversion");}

  static inline bool isPacketType() { 
    return SerializationTraits<int, char>::supportsDirectSerialization;
  }
  static inline bool equal(const char  a, const char  b) { return (a==b); }
  static inline bool lessThan(const char  a, const char  b) { return (a<b); }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<char, char> globalMinMax(
    char min, char max, const Comm<int> &comm) { 
       return z2GlobalMinMax(min, max, comm);}
  static std::pair<char, char> minMax(const char *values, size_t n) { 
    return z2LocalMinMax(values, n);}
  static bool areConsecutive(const char *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<short> {
  typedef unsigned long ulong;
  static inline int hashCode(const short  a) { return static_cast<int>(a);}
  static inline double key(const short a){ 
    return static_cast<double>(a); }
  static short keyToGid(const double key){ 
    return static_cast<short>(key);}
  static inline std::string name()   { return("short");}
  static std::string stringify(short val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const short x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const short x, short &y) {y = x;}
  static inline void castTo(const short x, int &y) {y = int(x);}
  static inline void castTo(const short x, unsigned &y) {y = unsigned(x);}
  static inline void castTo(const short x, long &y) {y = long(x);}
  static inline void castTo(const short x, ulong &y) {y = ulong(x);}
#ifdef HAVE_LONG_LONG
  typedef long long llong;
  typedef unsigned long long ullong;
  static inline void castTo(const short x, llong &y) {y = llong(x);}
  static inline void castTo(const short x, ullong &y) {y = ullong(x);}
#endif
  template <typename T1, typename T2>
    static inline void castTo(const short x, std::pair<T1, T2> &y) { 
      throw std::runtime_error("invalid conversion"); }

  static inline bool isPacketType() { 
    return SerializationTraits<int, short>::supportsDirectSerialization;
  }
  static inline bool equal(const short a, const short b) { 
    return (a==b) ; }
  static inline bool lessThan(const short a, const short b) { 
    return (a<b) ; }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<short, short> globalMinMax(
    short min, short max, const Comm<int> &comm) { 
      return z2GlobalMinMax(min, max, comm);}
  static std::pair<short, short> minMax(
    const short *values, size_t n) { return z2LocalMinMax(values, n);}
  static bool areConsecutive(const short *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<int> {
  typedef unsigned long ulong;
  static inline int hashCode(const int a) { return a; }
  static inline double key(const int a){ return static_cast<double>(a); }
  static int keyToGid(const double key){ return static_cast<int>(key); }
  static inline std::string name()    { return("int");}
  static std::string stringify(int val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return true; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const int x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const int x, short &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const int x, int &y) {y = x;}
  static inline void castTo(const int x, unsigned &y) {y = unsigned(x);}
  static inline void castTo(const int x, long &y) {y = long(x);}
  static inline void castTo(const int x, ulong &y) {y = ulong(x);}
#ifdef HAVE_LONG_LONG
  typedef long long llong;
  typedef unsigned long long ullong;
  static inline void castTo(const int x, llong &y) {y = llong(x);}
  static inline void castTo(const int x, ullong &y) {y = ullong(x);}
#endif
  template <typename T1, typename T2>
    static inline void castTo(const int x, std::pair<T1, T2>  &y) { 
      throw std::runtime_error("invalid conversion"); }

  static inline bool isPacketType() { 
    return SerializationTraits<int, int>::supportsDirectSerialization;
  }
  static inline bool equal(const  int a, const  int b) { 
    return (a==b) ; }
  static inline bool lessThan(const  int a, const  int b) { 
    return (a<b) ; }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<int, int> globalMinMax(
    int min, int max, const Comm<int> &comm) { 
    return z2GlobalMinMax(min, max, comm);}
  static std::pair<int, int> minMax(const int *values, size_t n) { 
    return z2LocalMinMax(values, n);}
  static bool areConsecutive(const int *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned> {
  typedef unsigned long ulong;
  static inline int hashCode(const unsigned a) { 
    return static_cast<int>(a); }
  static inline double key(const unsigned a){ 
    return static_cast<double>(a); }
  static unsigned keyToGid(const double key){  
    return static_cast<unsigned>(key);}
  static inline std::string name()             { return("unsigned");}
  static std::string stringify(unsigned val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const unsigned x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const unsigned x, short &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const unsigned x, int &y) {y = int(x);}
  static inline void castTo(const unsigned x, unsigned &y) {y = x;}
  static inline void castTo(const unsigned x, long &y) {y = long(x);}
  static inline void castTo(const unsigned x, ulong &y) {y = ulong(x);}
#ifdef HAVE_LONG_LONG
  typedef unsigned long long ullong;
  typedef long long llong;
  static inline void castTo(const unsigned x, llong &y) {y = llong(x);}
  static inline void castTo(const unsigned x, ullong &y) {y = ullong(x);}
#endif
  template <typename T1, typename T2>
    static inline void castTo(const unsigned x, std::pair<T1, T2>  &y) { 
      throw std::runtime_error("invalid conversion"); }
      
  static inline bool isPacketType() {
    return SerializationTraits<int, unsigned>::supportsDirectSerialization;
  }
  static inline bool equal(const unsigned a, const unsigned b) { 
    return (a==b) ; }
  static inline bool lessThan(const unsigned a, const unsigned b) { 
    return (a<b) ; }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<unsigned, unsigned> globalMinMax(
    unsigned min, unsigned max, const Comm<int> &comm) { 
      return z2GlobalMinMax(min, max, comm);}
  static std::pair<unsigned, unsigned> minMax(
    const unsigned *values, size_t n) { return z2LocalMinMax(values, n);}
  static bool areConsecutive(const unsigned *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<long> {
  typedef unsigned long ulong;
  static inline int hashCode(const long a) { 
    unsigned total=0;
    for (unsigned i=0, bits=0; i < sizeof(long); i++, bits += 8){
      total += static_cast<unsigned>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const long a){ return static_cast<double>(a); }
  static long keyToGid(const double key){ return static_cast<long>(key); }
  static inline std::string name()    { return("long");}
  static std::string stringify(long val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const long x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const long x, short &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const long x, int &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const long x, unsigned &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const long x, long &y) {y = x;}
  static inline void castTo(const long x, ulong &y) {y = ulong(x);}
#ifdef HAVE_LONG_LONG
  typedef unsigned long long ullong;
  typedef long long llong;
  static inline void castTo(const long x, llong &y) {y = llong(x);}
  static inline void castTo(const long x, ullong &y) {y = ullong(x);}
#endif
  template <typename T1, typename T2>
    static inline void castTo(const long x, std::pair<T1, T2>  &y) { 
      throw std::runtime_error("invalid conversion"); }
      
  static inline bool isPacketType() { 
    return SerializationTraits<int, long>::supportsDirectSerialization;
  }
  static inline bool equal(const long a, const long b) { 
    return (a==b) ; }
  static inline bool lessThan(const long a, const long b) { 
    return (a<b) ; }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<long, long> globalMinMax(
    long min, long max, const Comm<int> &comm) { 
      return z2GlobalMinMax(min, max, comm);}
  static std::pair<long, long> minMax(const long *values, size_t n) { 
    return z2LocalMinMax(values, n);}
  static bool areConsecutive(const long *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned long> {
  typedef unsigned long ulong;
  static inline int hashCode(const ulong a) { 
    return IdentifierTraits<long>::hashCode(static_cast<ulong>(a)); }
  static inline double key(const ulong a){ 
    return static_cast<double>(a);}
  static ulong keyToGid(const double key){ 
    return static_cast<ulong>(key);}
  static inline std::string name()   { return("unsigned long");}
  static std::string stringify(ulong val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const ulong x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ulong x, short &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ulong x, int &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ulong x, unsigned &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ulong x, long &y) {y = long(x);}
  static inline void castTo(const ulong x, ulong &y) {y = x;}
#ifdef HAVE_LONG_LONG
  typedef unsigned long long ullong;
  typedef long long llong;
  static inline void castTo(const ulong x, llong &y) {y = llong(x);}
  static inline void castTo(const ulong x, ullong &y) {y = ullong(x);}
#endif
  template <typename T1, typename T2>
    static inline void castTo(const ulong x, std::pair<T1, T2>  &y) { 
      throw std::runtime_error("invalid conversion"); }
      
  static inline bool isPacketType() { 
    return 
      SerializationTraits<int, ulong>::supportsDirectSerialization;
  }
  static inline bool equal( const unsigned long a, 
    const unsigned long b) { return (a==b) ; }
  static inline bool lessThan( const unsigned long a, 
    const unsigned long b) { return (a<b) ; }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<unsigned long, unsigned long> globalMinMax(
    unsigned long min, unsigned long max, const Comm<int> &comm) { 
      return z2GlobalMinMax(min, max, comm);}
  static std::pair<unsigned long, unsigned long> minMax(const 
    unsigned long *values, size_t n) {return z2LocalMinMax(values, n);}
  static bool areConsecutive(const unsigned long *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

#ifdef HAVE_LONG_LONG

template<>
struct IdentifierTraits<long long> {
  typedef unsigned long ulong;
  typedef unsigned long long ullong;
  typedef long long llong;
  static inline int hashCode(const llong a) { 
    unsigned total=0;
    for (unsigned i=0, bits=0; i < sizeof(llong); i++, bits += 8){
      total += static_cast<unsigned>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const llong a){ 
    return static_cast<double>(a);}
  static llong keyToGid(const double key){ 
    return static_cast<llong>(key);}
  static inline std::string name()    { return("llong");}
  static std::string stringify(llong val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return true; }

  static inline void castTo(const llong x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const llong x, short &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const llong x, int &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const llong x, unsigned &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const llong x, long &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const llong x, ulong &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const llong x, llong &y) {y = x;}
  static inline void castTo(const llong x, ullong &y) {y = ullong(x);}
  template <typename T1, typename T2>
    static inline void castTo(const llong x, std::pair<T1, T2>  &y) { 
      throw std::runtime_error("invalid conversion"); }
      
  static inline bool isPacketType() { 
     return SerializationTraits<int, llong>::supportsDirectSerialization;
  }
  static inline bool equal( const llong a, const llong b) { 
    return (a==b) ; }
  static inline bool lessThan( const llong a, const llong b) { 
    return (a<b) ; }
  static inline bool is_valid_id_type() { return true; }
  static std::pair<llong, llong> globalMinMax(
    llong min, llong max, const Comm<int> &comm) { 
      return z2GlobalMinMax(min, max, comm);}
  static std::pair<llong, llong> minMax(
    const llong *values, size_t n) { return z2LocalMinMax(values, n);}
  static bool areConsecutive(const llong *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned long long > {
  typedef unsigned long ulong;
  typedef unsigned long long ullong;
  typedef long long llong;
  static inline int hashCode(const ullong a) { 
    return IdentifierTraits<long long>::hashCode(static_cast<long long>(a));}
  static inline double key(const ullong a){ 
    return static_cast<double>(a); }
  static ullong keyToGid(const double key){ 
    return static_cast<ullong>(key);}
  static inline std::string name()    { return("unsigned long long");}
  static std::string stringify(ullong val) { return stringifyOrdinal(val);}
  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return false; }

  static inline void castTo(const ullong x, char &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ullong x, short &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ullong x, int &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ullong x, unsigned &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ullong x, long &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ullong x, ulong &y) {
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const ullong x, llong &y) {y = llong(x);}
  static inline void castTo(const ullong x, ullong &y) {y = x;}
  template <typename T1, typename T2>
    static inline void castTo(const ullong x, std::pair<T1, T2>  &y) { 
      throw std::runtime_error("invalid conversion"); }
      
  static inline bool isPacketType() { 
    return 
      SerializationTraits<int,ullong>::supportsDirectSerialization;
  }
  static inline bool equal(const ullong a, const ullong b) {return (a==b);}
  static inline bool lessThan(const ullong a, const ullong b) {return (a<b);}
  static inline bool is_valid_id_type() { return true; }
  static std::pair<ullong, ullong> globalMinMax(
    ullong min, ullong max, const Comm<int> &comm) { 
      return z2GlobalMinMax(min, max, comm);}
  static std::pair<ullong, ullong> minMax(const 
    ullong *values, size_t n) { return z2LocalMinMax(values, n);}
  static bool areConsecutive(const ullong *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

#endif

template<typename T1, typename T2>
struct IdentifierTraits<std::pair<T1, T2> > {
  typedef unsigned long ulong;
  typedef unsigned long long ullong;
  typedef long long llong;
  typedef std::pair<T1, T2> pair_t;
  typedef typename std::pair<pair_t, pair_t> pairPair_t;

  static inline int hashCode(const pair_t p)  {
    return IdentifierTraits<T1>::hashCode(p.first) + 
      IdentifierTraits<T2>::hashCode(p.second);
  }

  static inline double key(const pair_t p)  {
    int nbits = sizeof(T1)*8;
    intmax_t part1(p.first);
    intmax_t part2(p.second);
    return static_cast<double>((part2 << nbits) | part1);
  }

  static pair_t keyToGid(const double key){
    intmax_t t1_size_mask = 0xff;
    for (unsigned n=1; n < sizeof(T1); n++){
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

  static std::string stringify(pair_t val) { 
    std::ostringstream oss;
    oss << "pair<" << IdentifierTraits<T1>::name();
    oss << "," << IdentifierTraits<T2>::name();
    oss << ">(" << val.first << "," << val.second << ")";
    return oss.str();
  }

  static inline bool isHashKeyType() { return false; }
  static inline bool isGlobalOrdinal() { return false; }

  static inline void castTo(const pair_t x, char &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, short &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, int &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, unsigned &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, long &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, ulong &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, llong &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, ullong &y) { 
      throw std::runtime_error("invalid conversion"); }
  static inline void castTo(const pair_t x, pair_t &y) { y = x;}

  static inline bool isPacketType() { 
    return 
      SerializationTraits<int,pair_t >::supportsDirectSerialization;
  }
  static inline bool equal( const pair_t a, 
    const pair_t b) { 
    return ((a.first==b.first) && (a.second==b.second)); }
  static inline bool lessThan( const pair_t a, 
    const pair_t b) { 
      throw std::logic_error("invalid call");
      return false;}
  static inline bool is_valid_id_type() { 
    return (sizeof(T1)+sizeof(T2) <= sizeof(intmax_t)); }
  static pairPair_t globalMinMax(
    pair_t min, pair_t max, const Comm<int> &comm) { 
      throw std::logic_error("invalid call");
      return pairPair_t(); }
  static pairPair_t minMax(const pair_t *values, size_t n) { 
      throw std::logic_error("invalid call");
      return pairPair_t(); }
  static bool areConsecutive(const pair_t *val, size_t n){ 
      throw std::logic_error("invalid call");
      return false; }
};

/*! \endcond
 */

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
