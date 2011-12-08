// @HEADER
// ***********************************************************************
//            copyright
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERTRAITS
#define _ZOLTAN2_IDENTIFIERTRAITS

#include <Zoltan2_Standards.hpp>

#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_HashUtils.hpp>

#include <utility>
#include <iostream>
#include <sstream>
#include <string>

using Teuchos::SerializationTraits;

/*! \file Zoltan2_IdentifierTraits.hpp
  \brief Defines basic traits for application supplied global IDs.

  The data types permitted for global identifiers for Zoltan2 callers include
  some that are not represented in Teuchos::OrdinalTraits.  A common case is
  when a matrix nonzero is represented as an (i,j) pair.

  Zoltan2 uses the IdentifierTraits structures to manage the application
  supplied identifiers.
*/

/*! \namespace Zoltan2
  \brief Internal Zoltan2 namespace.

  This namespace contains all symbols in the Zoltan2 library.
*/
namespace Zoltan2
{


/*! \brief helper function to find min and max of values
 */
template <typename T>
  std::pair<T, T> z2LocalMinMax(const T *val, size_t n)
{
  if (n < 1) return std::pair<T,T>(0,0);  // TODO

  T min = val[0], max = val[0];
  for (size_t i=1; i < n; i++){
    if (val[i] < min) min = val[i];
    else if (val[i] > max) max = val[i];
  }
  return std::pair<T,T>(min,max);
}

/*! \brief helper function to determine if values are consecutive
 */
template <typename T>
  bool z2AreConsecutive(const T *val, size_t n)
{
  if (n == 0) return true;

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

/*! \brief helper function write value to a string
 */
template <typename T>
  std::string stringifyOrdinal(T ordinal)
{
  std::ostringstream oss;
  oss << ordinal;
  return oss.str();
}

/*!  \brief Structure to catch invalid Indentifier types.
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
      Assumption: Any integer value can fit in a double.
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

  /*! \brief Determine if two identifiers of the same type are equal.
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

  /*! \brief Compute b - a, if possible
      \result b-a if a and b are ordinals,
                 throw an error otherwise.
   */
  static inline T difference(const T a, const T b) {
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

/*! \cond IdentifierTraits_definitions
 */

template<>
struct IdentifierTraits<char> {
  typedef char T;
  static inline int hashCode(const T c) {return static_cast<int>(c);}
  static inline double key(const T c){return static_cast<double>(c); }
  static T keyToGid(const double key){return static_cast<T>(key); }
  static inline std::string name()     {return("char");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() { 
   return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T  a, const T  b) {return (a==b); }
  static inline bool lessThan(const T  a, const T  b) {return (a<b); }
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(const T *values, size_t n) { 
   return z2LocalMinMax(values, n);}
  static bool areConsecutive(const T *val, size_t n){ 
   return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<int16_t> {
  typedef int16_t T;
  static inline int hashCode(const T  a) {return static_cast<int>(a);}
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("int16_t");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() { 
  return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T a, const T b) {return (a==b);}
  static inline bool lessThan(const T a, const T b) {return (a<b);}
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(
    const T *values, size_t n) {return z2LocalMinMax(values, n);}
  static bool areConsecutive(const T *val, size_t n){ 
  return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<uint16_t> {
  typedef uint16_t T;
  static inline int hashCode(const T  a) {return static_cast<int>(a);}
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("uint16_t");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() { 
  return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T a, const T b) {return (a==b);}
  static inline bool lessThan(const T a, const T b) {return (a<b);}
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(
    const T *values, size_t n) {return z2LocalMinMax(values, n);}
  static bool areConsecutive(const T *val, size_t n){ 
  return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<int32_t> {
  typedef int32_t T;
  static inline int hashCode(const T a) {return static_cast<int>(a); }
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key); }
  static inline std::string name()    {return("int32_t");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return true; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() { 
    return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T a, const  T b) {return (a==b);}
  static inline bool lessThan(const T a, const T b) {return (a<b);}
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(const T *values, size_t n) { 
    return z2LocalMinMax(values, n);}
  static bool areConsecutive(const T *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<uint32_t> {
  typedef uint32_t T;
  static inline int hashCode(const T a) {return static_cast<int>(a); }
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()             {return("uint32_t");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() {
    return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T a, const T b) {return (a==b);}
  static inline bool lessThan(const T a, const T b) {return (a<b);}
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(const T *vals, size_t n){ 
    return z2LocalMinMax(vals, n);}
  static bool areConsecutive(const T *val, size_t n){ 
   return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<int64_t> {
  typedef int64_t T;
  static inline int hashCode(const T a) { 
    unsigned total=0;
    for (unsigned i=0, bits=0; i < sizeof(T); i++, bits += 8){
      total += static_cast<unsigned>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key); }
  static inline std::string name()    {return("int64_t");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() { 
    return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T a, const T b) {return (a==b);}
  static inline bool lessThan(const T a, const T b) {return (a<b);}
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(const T *values, size_t n) { 
    return z2LocalMinMax(values, n);}
  static bool areConsecutive(const T *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<uint64_t> {
  typedef uint64_t T;
  static inline int hashCode(const T a) { 
    return IdentifierTraits<int64_t>::hashCode(static_cast<int64_t>(a)); }
  static inline double key(const T a){return static_cast<double>(a);}
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("uint64_t");}
  static std::string stringify(T val) {return stringifyOrdinal(val);}
  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return true; }
  static inline bool isPacketType() { 
    return SerializationTraits<int, T>::supportsDirectSerialization;
  }
  static inline bool equal(const T a, const T b) {return (a==b);}
  static inline bool lessThan(const T a, const T b) {return (a<b);}
  static inline T difference(const T a, const T b) { return (b-a); }
  static inline bool is_valid_id_type() {return true; }
  static std::pair<T, T> minMax(const T *values, size_t n)
    {return z2LocalMinMax(values, n);}
  static bool areConsecutive(const T *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<typename T1, typename T2>
struct IdentifierTraits<std::pair<T1, T2> > {
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

  static inline bool isHashKeyType() {return false; }
  static inline bool isGlobalOrdinal() {return false; }

  static inline bool isPacketType() { 
    return SerializationTraits<int,pair_t >::supportsDirectSerialization;
  }

  static inline bool equal( const pair_t a, const pair_t b) { 
    return ((a.first==b.first) && (a.second==b.second)); }

  static inline bool lessThan( const pair_t a, const pair_t b) { 
      throw std::logic_error("invalid call");
      return false;}

  static inline pair_t difference( const pair_t a, const pair_t b) { 
      throw std::logic_error("invalid call");
      return false;}

  static inline bool is_valid_id_type() { 
    return (sizeof(T1)+sizeof(T2) <= sizeof(int64_t)); }

  static pairPair_t minMax(const pair_t *values, size_t n) { 
      throw std::logic_error("invalid call");
      return pairPair_t(); }

  static bool areConsecutive(const pair_t *val, size_t n){ 
      throw std::logic_error("invalid call");
      return false; }
};

//////////////////////////////////////////////////////////////
//  A helper function.  Are the ordinals globally consecutive?
//    If IDs are globally consecutive, then return true and
//    set nprocs+1 length dist array to proc starting IDs.
//////////////////////////////////////////////////////////////

template <typename T>
  bool globallyConsecutiveOrdinals(const T* val, size_t len, size_t globalLen, 
   const Comm<int> &comm, const Environment &env, ArrayView<T> dist)
{
  bool globallyConsecutive = false;

  if (!IdentifierTraits<T>::isGlobalOrdinal()){
    return globallyConsecutive;
  }

  // Are they consecutive and increasing with process rank?

  bool locallyConsecutiveIncreasing =
    IdentifierTraits<T>::areConsecutive(val, len);

  int localFlag = (locallyConsecutiveIncreasing ? 1 : 0);
  int globalFlag = 0;

  try{
    reduceAll<int, int>(comm, Teuchos::REDUCE_MIN, 1, &localFlag, &globalFlag);
  }
  catch (const std::exception &e) {
    Z2_THROW_OUTSIDE_ERROR(env, e);
  }

  if (globalFlag == 1){

    int64_t lMin = INT64_MAX, gMin;
    int64_t lMax = INT64_MIN, gMax;

    if (len > 0){
      lMin = Teuchos::as<int64_t>(val[0]);
      lMax = Teuchos::as<int64_t>(val[len-1]);
    }

    try{
      reduceAll<int, int64_t>(comm, Teuchos::REDUCE_MIN, 1, &lMin, &gMin);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(env, e);
    }

    try{
      reduceAll<int, int64_t>(comm, Teuchos::REDUCE_MAX, 1, &lMax, &gMax);
    }
    catch (const std::exception &e) {
      Z2_THROW_OUTSIDE_ERROR(env, e);
    }

    if (gMax - gMin + 1 == globalLen){

      int nprocs = comm.getSize();

      Array<int64_t> sendBuf(1, lMin);
      Array<int64_t> recvBuf(nprocs);

      try{
        Teuchos::gatherAll<int, int64_t>(comm, 1, sendBuf.getRawPtr(), nprocs, 
          recvBuf.getRawPtr());
      }
      Z2_FORWARD_EXCEPTIONS;

      globallyConsecutive = true;

      for (int i=1; i < nprocs; i++){
        if (recvBuf[i-1] < recvBuf[i])
          continue;
        globallyConsecutive = false;
        break;
      }

      if (globallyConsecutive){
        if (dist.size() <= nprocs)
          throw std::logic_error("dist not preallocated");

        for (int i=0; i < nprocs; i++){
          dist[i] = static_cast<T>(recvBuf[i]);
        }
 
        dist[nprocs] = globalLen + dist[0];
      }
    }
  }
  return globallyConsecutive;
}

/*! \endcond
 */

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
