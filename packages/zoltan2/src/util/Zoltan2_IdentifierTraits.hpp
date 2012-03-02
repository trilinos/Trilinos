// @HEADER
// ***********************************************************************
//            copyright
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IdentifierTraits.hpp
   \brief Defines basic traits for user global identifiers.
*/

#ifndef _ZOLTAN2_IDENTIFIERTRAITS
#define _ZOLTAN2_IDENTIFIERTRAITS

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_AlltoAll.hpp>

#include <Teuchos_SerializationTraits.hpp>
#include <Teuchos_HashUtils.hpp>

#include <utility>
#include <iostream>
#include <sstream>
#include <string>

using Teuchos::SerializationTraits;

namespace Zoltan2
{

/*! \brief helper function to find min and max of array of user Ids
 */
template <typename T>
  std::pair<T, T> z2LocalMinMax(const T *val, size_t n)
{
  if (n < 1) return std::pair<T,T>(0,0);

  T min = val[0], max = val[0];
  for (size_t i=1; i < n; i++){
    if (val[i] < min) min = val[i];
    else if (val[i] > max) max = val[i];
  }
  return std::pair<T,T>(min,max);
}

/*! \brief helper function to find global min and max of array of user Ids
 */
template <typename T>
  void z2GlobalMinMax(const Comm<int> &comm, 
    const T &localMin, const T &localMax, T &globalMin, T &globalMax)
{
  int nprocs = comm.getSize();

  if (nprocs < 2){
    globalMin = localMin;
    globalMax = localMax;
    return;
  }
  
  // Need an environment for error messages
  ParameterList params;
  RCP<Comm<int> > comm2 = comm.duplicate();
  RCP<const Comm<int> > comm2Const = rcp_const_cast<const Comm<int> >(comm2);
  Environment env(params, comm2Const);

  ArrayRCP<T> recvBufMin;
  ArrayRCP<T> recvBufMax;
  Array<T> sendBufMin(nprocs, localMin);
  Array<T> sendBufMax(nprocs, localMax);

  // Must use Zoltan2::AlltoAll because T may not be a valid
  //   Teuchos packet type.

  try{
    AlltoAll<T,int>(comm, env, sendBufMin.view(0,nprocs), 1, recvBufMin);
  }
  Z2_FORWARD_EXCEPTIONS; 

  try{
    AlltoAll<T,int>(comm, env, sendBufMax.view(0,nprocs), 1, recvBufMax);
  }
  Z2_FORWARD_EXCEPTIONS; 

  globalMin = recvBufMin[0];
  globalMax = recvBufMax[0];

  for (int i=0; i < nprocs; i++){
    if (recvBufMin[i] < globalMin) globalMin = recvBufMin[i];
    if (recvBufMax[i] > globalMax) globalMax = recvBufMax[i];
  }
}

/*! \brief helper function to determine if list of user Ids are consecutive
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

/*! \brief helper function write a user ID to a string
 */
template <typename T>
  std::string stringifyOrdinal(T ordinal)
{
  std::ostringstream oss;
  oss << ordinal;
  return oss.str();
}

/*!  \brief Structure to catch invalid user ID types.
 */
template<typename T>
struct UndefIdTraits
{
static inline T invalid() { 
  return T::UsingInvalidGlobalIdentifierDataType(); 
}
};


/*!  \brief The functions to be defined for users' global ID types.

  The template parameter is the user's global ID data type.

  The data types permitted for global identifiers for Zoltan2 callers 
  may include those that are not represented in Teuchos::OrdinalTraits.  
  A common case is
  when a matrix nonzero is represented as an (i,j) pair.

  In such a case, Zoltan2 will map them to
  a list of new IDs that \em are Teuchos Ordinals.  All computation will
  be in the space of the new global numbers.  When the Solution object
  is written, the internal global numbers are mapped back to the user's
  global IDs.

  During this process of determining if the user's IDs can be used by
  Zoltan2, and mapping them to new global numbers if they can not,
  the Traits defined here are used to manipulate the user's IDs.

  Traits are defined for the following types:
    \li char
    \li unsigned char
    \li short
    \li unsigned short
    \li int
    \li unsigned int
    \li long
    \li unsigned long
    \li long long
    \li unsigned long long
    \li std::pair<T1, T2>

  The <tt> long long </tt>  and <tt> unsigned long long </tt> traits are only 
  defined if Trilinos was configured with the TEUCHOS_ENABLE_LONG_LONG_INT
  option.

  The translation mechanism relies on the ability to assign each
  global ID a locally unique key, which is type \c double.  For this reason
  we assume that if the global ID is integral, its size is not greater
  than the size of a double.  And, in the case of \c std::pair<T1, T2>,
  the sum of the sizes of \c T1 and \c T2 is not greater than the
  size of a double.

  If the user's global ID type does not appear in the above list, it
  can be added by the user in his or her application code.  See the example
  in tobeWritten.cpp.

  \todo write an example where user's global ID is a C-struct containing
        \c i and \c j indices.

  Developer note: By convention we use \c gid_t as the users global ID
  data type and \c gno_t as the data type used internally by Zoltan2.
*/

template<typename T>
struct IdentifierTraits {

  /*! \brief Compute an integer hash code for the user's global ID.

      \param id  the user's global id
      \result the integer code, which need not be unique for each id
   */
  static int hashCode(const T id) {
   return UndefIdTraits<int>::invalid();
  }

  /*! \brief Compute a key which will be unique for each id.

      \param id  the user's global id
      \result the key
   */
  static double key(const T id){
   return UndefIdTraits<double>::invalid();
  }

  /*! \brief Convert a key back to the identifier that generated it.

      \param x the key
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
		
  /*! \brief A string displaying the value.

      \param  val  the value to represent as a string
      \result The string
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

     Data types are those with a definition in Teuchos::OrdinalTraits.
   */
  static inline bool isGlobalOrdinal() {
    return UndefIdTraits<bool>::invalid();
  }

  /*!  \brief Determine whether the data type can be used in 
                 Teuchos communication.

      \result true if it can be a Teuchos Packet type

  Packet data types are those with a definition in Teuchos::SerializationTraits.
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

      \param a  The \em a of a<b
      \param b  The \em b of a<b
      \result true if a < b, false if !(a<b)

       A \c std::logic_error is throw at runtime if T is a type
            that can not be ordered.
   */
  static inline bool lessThan(const T a, const T b) {
    return UndefIdTraits<bool>::invalid();
  }

  /*! \brief Compute b - a, if possible

      \param a  The \em a of b-a
      \param b  The \em b of b-a
      \result the value b-a 

       A \c std::logic_error is throw at runtime if the operation
                 is not valid for T.
   */
  static inline T difference(const T a, const T b) {
    return UndefIdTraits<bool>::invalid();
  }

  /*! \brief Determine if the data type is one for which IdentifierTraits are 
                   defined

      \result true if data type has definition
   */
  static inline bool is_valid_id_type() { return false; }

  /*! \brief Return the minimum and maximum of a list of values.
      \param values   a pointer to \c numValues values
      \param numValues   the number of values to consider
      \result A pair with the minimum value followed by the maximum value 

       A \c std::logic_error is throw at runtime if T is a type
            that can not be ordered.
   */
  static std::pair<T, T> minMax(const T *values, size_t numValues) {
    return UndefIdTraits<std::pair<T, T> >::invalid();
  }

  /*! \brief Find global minimum and maximum

      \param comm  Communicator for global operation
      \param localMin  The local minimum
      \param localMax  The local maximum
      \param globalMin  On return, the global minimum
      \param globalMax  On return, the global maximum

       A \c std::logic_error is throw at runtime if T is a type
            that can not be ordered.
   */
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax)
        {UndefIdTraits<std::pair<T, T> >::invalid();}

  /*! \brief Determine if the values are locally increasing consecutive
      
      \param val  a pointer to \c n values
      \param n  the number of values in the list
      \result true if the values are increasing consecutive,
            false if otherwise

       A \c std::logic_error is throw at runtime if T is a type
            that can not be ordered.
   */
  static bool areConsecutive(const T *val, size_t n){
    return UndefIdTraits<bool>::invalid();
  }

};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
   return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned char> {
  typedef unsigned char T;
  static inline int hashCode(const T c) {return static_cast<int>(c);}
  static inline double key(const T c){return static_cast<double>(c); }
  static T keyToGid(const double key){return static_cast<T>(key); }
  static inline std::string name()     {return("unsigned char");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
   return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<short> {
  typedef short T;
  static inline int hashCode(const T  a) {return static_cast<int>(a);}
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("short");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
  return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned short> {
  typedef unsigned short T;
  static inline int hashCode(const T  a) {return static_cast<int>(a);}
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("unsigned short");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
  return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<int> {
  typedef int T;
  static inline int hashCode(const T  a) {return static_cast<int>(a);}
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("int");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
  return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned int> {
  typedef unsigned int T;
  static inline int hashCode(const T  a) {return static_cast<int>(a);}
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("unsigned int");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
  return z2AreConsecutive(val, n); }
};


template<>
struct IdentifierTraits<long> {
  typedef long T;
  static inline int hashCode(const T a) {
    unsigned total=0;
    for (unsigned i=0, bits=0; i < sizeof(T); i++, bits += 8){
      total += static_cast<unsigned>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key); }
  static inline std::string name()    {return("long");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned long> {
  typedef unsigned long T;
  static inline int hashCode(const T a) {
    return IdentifierTraits<long>::hashCode(static_cast<long>(a)); }
  static inline double key(const T a){return static_cast<double>(a);}
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("unsigned long");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){
    return z2AreConsecutive(val, n); }
};

#ifdef HAVE_ZOLTAN2_LONG_LONG

template<>
struct IdentifierTraits<long long> {
  typedef long long T;
  static inline int hashCode(const T a) {
    unsigned total=0;
    for (unsigned i=0, bits=0; i < sizeof(T); i++, bits += 8){
      total += static_cast<unsigned>((a & (0xff << bits) ) >> bits);
    }
    return static_cast<int>(total);
  }
  static inline double key(const T a){return static_cast<double>(a); }
  static T keyToGid(const double key){return static_cast<T>(key); }
  static inline std::string name()    {return("long long");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

template<>
struct IdentifierTraits<unsigned long long> {
  typedef unsigned long long T;
  static inline int hashCode(const T a) {
    return IdentifierTraits<long long>::hashCode(static_cast<long long>(a)); }
  static inline double key(const T a){return static_cast<double>(a);}
  static T keyToGid(const double key){return static_cast<T>(key);}
  static inline std::string name()   {return("unsigned long long");}
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
  static void globalMinMax(const Comm<int> &comm,
      const T &localMin, const T &localMax, T &globalMin, T &globalMax){
    z2GlobalMinMax(comm, localMin, localMax, globalMin, globalMax);}
  static bool areConsecutive(const T *val, size_t n){ 
    return z2AreConsecutive(val, n); }
};

#endif

template<typename T1, typename T2>
struct IdentifierTraits<std::pair<T1, T2> > {
  typedef std::pair<T1, T2> pair_t;
  typedef typename std::pair<pair_t, pair_t> pairPair_t;

  static inline int hashCode(const pair_t p)  {
    return IdentifierTraits<T1>::hashCode(p.first) +
      IdentifierTraits<T2>::hashCode(p.second);
  }

  static inline double key(const pair_t p)  {
    if (sizeof(T1) + sizeof(T2) > sizeof(double))
      throw std::runtime_error("pair gid is invalid");

    double keyVal;
    char *cx = reinterpret_cast<char *>(&keyVal);
    T1 *xpos = reinterpret_cast<T1 *>(cx);
    T2 *ypos = reinterpret_cast<T2 *>(cx + sizeof(T1));
    *xpos = p.first;
    *ypos = p.second;
    
    return keyVal;
  }

  static pair_t keyToGid(const double key){
    if (sizeof(T1) + sizeof(T2) > sizeof(double))
      throw std::runtime_error("pair gid is invalid");

    const char *cx = reinterpret_cast<const char *>(&key);
    const T1 *xpos = reinterpret_cast<const T1 *>(cx);
    const T2 *ypos = reinterpret_cast<const T2 *>(cx + sizeof(T1));

    return std::pair<T1,T2>(*xpos, *ypos);
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
    return (sizeof(T1)+sizeof(T2) <= sizeof(double)); }

  static pairPair_t minMax(const pair_t *values, size_t n) {
      throw std::logic_error("invalid call");
      return pairPair_t(); }

  static void globalMinMax(const Comm<int> &comm,
      const pair_t &localMin, const pair_t &localMax, 
      pair_t &globalMin, pair_t &globalMax){
      throw std::logic_error("invalid call");}

  static bool areConsecutive(const pair_t *val, size_t n){ 
      throw std::logic_error("invalid call");
      return false; }
};

//////////////////////////////////////////////////////////////
//  A helper function.  If T is an ordinal, are the values
//    globally consecutive?  If so return true, otherwise
//    return false.
//
//  On return, globalLen is set to the sum of the local lengths.
//
//  If T is an ordinal, but the list is not globally consecutive, 
//    on return dist[0] is set to the global minimum of
//    the values and dist[1] to the global maximum.
//    
//  If T is an ordinal and the list is globally consecutive,
//    on return dist[p] is set to val[0] on process p.  dist[nprocs]
//    is set to one past the global maximum value.
//////////////////////////////////////////////////////////////

template <typename T>
  bool globallyConsecutiveOrdinals(
    const Comm<int> &comm, const Environment &env, const T* val, size_t len,
    ArrayRCP<T> &dist, size_t &globalLen)
{
  try{
    reduceAll<int, size_t>(comm, Teuchos::REDUCE_SUM, 1, &len, &globalLen);
  }
  Z2_THROW_OUTSIDE_ERROR(env, e);

  if (!IdentifierTraits<T>::isGlobalOrdinal()){
    return false;
  }

  // Get global minimum and maximum

  T gMin, gMax;
  T v0 = val[0];
  T v1 = val[len-1];

  try{
    IdentifierTraits<T>::globalMinMax(comm, v0, v1, gMin, gMax);
  }
  Z2_FORWARD_EXCEPTIONS; 

  T *minMax = new T [2];
  minMax[0] = gMin;
  minMax[1] = gMax;
  dist = arcp<T>(minMax, 0, 2);

  size_t g0 = Teuchos::as<size_t>(gMin);
  size_t g1 = Teuchos::as<size_t>(gMax);
  bool globallyConsecutive = false;

  if (g1 - g0 + 1 == globalLen){

    size_t sentinel = g1 + 1;
    int nprocs = comm.getSize();
    bool locallyConsecutive = IdentifierTraits<T>::areConsecutive(val, len);

    if (locallyConsecutive && nprocs==1){
      dist[nprocs] = Teuchos::as<T>(sentinel);
      return true;
    }

    int lFlag = (locallyConsecutive ? 1 : 0);
    int gFlag = 0;

    try{
      reduceAll<int, int>(comm, Teuchos::REDUCE_MIN, 1, &lFlag, &gFlag);
    }
    Z2_THROW_OUTSIDE_ERROR(env, e);

    if (gFlag == 0)  // not all processes have consecutive values
      return false;

    Array<size_t> sendBuf(nprocs);
    ArrayRCP<size_t> recvBuf;

    for (int i=0; i < nprocs; i++)
      sendBuf[i] = Teuchos::as<size_t>(v0);

    try{
      AlltoAll<size_t, int>(comm, env, sendBuf, 1, recvBuf);
    }
    Z2_FORWARD_EXCEPTIONS;

    globallyConsecutive = true;
    for (int i=1; globallyConsecutive && i < nprocs; i++)
      if (recvBuf[i] < recvBuf[i-1]) globallyConsecutive = false;

    if (globallyConsecutive){
      T *idDist = new T [nprocs+1];
      for (int i=0; i < nprocs; i++)
        idDist[i] = Teuchos::as<T>(recvBuf[i]);
      idDist[nprocs] = Teuchos::as<T>(sentinel);
      dist = arcp(idDist, 0, nprocs+1);
    }
  }

  return globallyConsecutive;
}

#endif  // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
