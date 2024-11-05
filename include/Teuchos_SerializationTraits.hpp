// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Teuchos_SerializationTraits.hpp
/// \brief Teuchos::SerializationTraits and Teuchos::DirectSerializationTraits definitions.
///
#ifndef TEUCHOS_SERIALIZATION_TRAITS_HPP
#define TEUCHOS_SERIALIZATION_TRAITS_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include <climits> // SIZE_MAX, ULONG_MAX, etc.

#ifdef HAVE_TEUCHOSCORE_QUADMATH
#include <quadmath.h>
#endif // HAVE_TEUCHOSCORE_QUADMATH

#ifdef HAVE_TEUCHOS_QD
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#endif

namespace Teuchos {

/// \class UndefinedSerializationTraits
/// \brief Report an error if a specialization of \c SerializationTraits is missing.
///
/// This class reports a compile-time error if you attempt to
/// instantiate it.  We use this class to make it easy to detect a
/// missing specialization of \c SerializationTraits.
template<typename T>
struct UndefinedSerializationTraits {
  //! This function should not compile if there is an attempt to instantiate!
  static inline T notDefined() {return(T::this_type_is_missing_a_specialization());}
};


/** \class SerializationTraits
 * \brief Serialization traits class for types T that use value semantics.
 *
 * This traits class describes how to convert between arrays of T, and
 * arrays of char.  We call the process of converting from T to char
 * "serialization," and from char to T "deserialization."  Teuchos
 * uses serialization and deserialization mainly for implementing
 * distributed-memory message-passing communication in a generic way.
 *
 * \tparam Ordinal The same template parameter as that of \c Comm.
 *   The integer type used to count the number of packets sent and
 *   received.
 *
 * \tparam T The type of the objects that this class shows how to
 *   serialize.
 *
 * Teuchos defines specializations of this class for many commonly
 * used types in distributed-memory communication, such as char, int
 * (signed and unsigned), float, double, long double and std::pair<P1, P2> for
 * certain types P1 and P2.  Depending on your Trilinos build options,
 * other specializations may be defined as well, for example for long
 * long int, double-double and quad-double real types (dd_real
 * resp. qd_real), or certain std::complex<T> specializations.  If a
 * specialization of this class does not exist for your type T, you
 * may define your own specialization.
 *
 * \note Before defining specializations of this class, make sure that
 *   they do not duplicate specializations already present in
 *   PyTrilinos (see packages/PyTrilinos/src/Teuchos_Traits.i)
 *
 * There are two different serialization modes: direct and indirect.
 * "Direct" serialization means that you can convert directly between
 * an object of type T and an array of char, of a specific length
 * dependent only on the type T and not on the particular instance.
 * Specifically, it means you can
 *
 * 1. reinterpret_cast a pointer to an instance of T into an array of
 *    char (which array has length dependent only on the type T and
 *    not on the specific T instance),
 * 2. serialize the resulting array of char, and finally
 * 3. deserialize by reading in the array of char and doing a
 *    reinterpret_cast back into a T.
 *
 * "Indirect" serialization is defined as any serialization method
 * more general than that.  The \c supportsDirectSerialization class
 * Boolean tells you whether this specialization supports direct
 * serialization.
 *
 * SerializationTraits is used by classes such as \c
 * ValueTypeSerializationBuffer, \c ConstValueTypeSerializationBuffer,
 * \c ValueTypeDeserializationBuffer, and \c
 * ConstValueTypeDeserializationBuffer.
 */
template <typename Ordinal, typename T>
class SerializationTraits {
public:

  //! @name Serialization type selection
  //@{

  /// \brief Whether the type T supports direct serialization.
  ///
  /// See the class documentation for definitions of "direct" and
  /// "indirect" serialization.
  static const bool supportsDirectSerialization = false;

  //@}

  //! @name Direct serialization functions (not defined if supportsDirectSerialization==false)
  //@{

  /** \brief Return the number of bytes for <tt>count</tt> objects. */
  static Ordinal fromCountToDirectBytes(const Ordinal count) {
    (void)count;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Convert the pointer type to <tt>char*</tt>. */
  static char* convertToCharPtr( T* ptr ) {
    (void)ptr;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Convert the pointer type to <tt>const char*</tt>. */
  static const char* convertToCharPtr( const T* ptr ) {
    (void)ptr;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
  static Ordinal fromDirectBytesToCount(const Ordinal bytes) {
    (void)bytes;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Convert the pointer type from <tt>char*</tt>. */
  static T* convertFromCharPtr( char* ptr ) {
    (void)ptr;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Convert the pointer type from <tt>char*</tt>. */
  static const T* convertFromCharPtr( const char* ptr ) {
    (void)ptr;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  //@}

  //! @name Indirect serialization functions (always defined and supported)
  //@{

  /** \brief Return the number of bytes for <tt>count</tt> objects. */
  static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                          const T buffer[]) {
    (void)count; (void)buffer;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Serialize to an indirect <tt>char[]</tt> buffer.
   *
   * \param  count
   *           [in] The number of objects to serialize.
   * \param  buffer
   *           [in] The objects to serialize.
   * \param  bytes
   *           [in] Number of bytes in <tt>charBuffer[]</tt>
   * \param  charBuffer
   *           [out] Array (length <tt>bytes</tt>) containing the serialized objects.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>bytes==fromCountToIndirectBytes(count)</tt>
   * </ul>
   */
  static void serialize (const Ordinal count,
                         const T buffer[],
                         const Ordinal bytes,
                         char charBuffer[])
  {
    (void)count; (void)buffer; (void)bytes; (void)charBuffer;
    UndefinedSerializationTraits<T>::notDefined();
  }

  /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
  static Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                          const char charBuffer[]) {
    (void)bytes; (void)charBuffer;
    UndefinedSerializationTraits<T>::notDefined();
    return 0;
  }

  /** \brief Deserialize from an indirect <tt>char[]</tt> buffer.
   *
   * \param  bytes
   *           [in] Number of bytes in <tt>charBuffer[]</tt>
   * \param  charBuffer
   *           [in] Array (length <tt>bytes</tt>) containing the serialized objects.
   * \param  count
   *           [in] The number of objects to deserialize.
   * \param  buffer
   *           [out] The deserialized objects.

   * <b>Preconditions:</b><ul>
   * <li><tt>count==fromIndirectBytesToCount(bytes)</tt>
   * </ul>
   */
  static void deserialize (const Ordinal bytes,
                           const char charBuffer[],
                           const Ordinal count,
                           T buffer[])
  {
    (void)bytes; (void)charBuffer; (void)count; (void)buffer;
    UndefinedSerializationTraits<T>::notDefined();
  }

  //@}

};

/** \class ValueTypeSerializer
 * \brief Serialization class for types T that use value semantics.
 *
 * This works similarly to SerializationTraits, except that it provides
 * a class that can be specialized for types T to implement serialization
 * through an object instead of a traits class.  Some types need other data
 * to help them serialize, and that data can be encapsulated into a non-static
 * specialization of this class.
 *
 * The default implementation is just given by the SerializationTraits for the
 * type T, and thus this class will be defined for any type that defines its
 * SerializationTraits.
 */
template <typename Ordinal, typename T>
class ValueTypeSerializer : public Teuchos::SerializationTraits<Ordinal,T> {};

/// \class DirectSerializationTraits
/// \brief Serialization traits for objects that support direct serialization.
///
/// "Direct" serialization means that you can convert directly
/// between an object of type T and an array of char, of a specific
/// length dependent only on the type T and not on the particular
/// instance.  Specifically, it means you can
///
/// 1. reinterpret_cast a pointer to an instance of T into an
///    array of char (which array has length dependent only on
///    the type T and not on the specific T instance),
/// 2. serialize the resulting array of char, and finally
/// 3. deserialize by reading in the array of char and doing a
///    reinterpret_cast back into a T.
///
/// "Indirect" serialization is defined as any serialization method
/// more general than that.
///
/// We use partial specializations of DirectSerializationTraits
/// (specialized on certain T types, not Ordinal) as public base
/// classes for the corresponding SerializationTraits specialization.
/// This provides high-performance default implementations of
/// serialization for commonly used types T (including char, int, and
/// double).
///
/// \tparam Ordinal The same template parameter as that of \c Comm.
///   The integer type used to count the number of packets sent and
///   received.
///
/// \tparam T The type of the objects that this class shows how to
///   serialize.
///
template <typename Ordinal, typename T>
class DirectSerializationTraits {
public:
  static const bool supportsDirectSerialization = true;
  // Direct serialization
  static Ordinal fromCountToDirectBytes(const Ordinal count)
    { return sizeof(T)*count; }
  static char* convertToCharPtr( T* ptr )
    { return reinterpret_cast<char*>(ptr); }
  static const char* convertToCharPtr( const T* ptr )
    { return reinterpret_cast<const char*>(ptr); }
  static Ordinal fromDirectBytesToCount(const Ordinal count)
    { return count/sizeof(T); }
  static T* convertFromCharPtr( char* ptr )
    { return reinterpret_cast<T*>(ptr); }
  static const T* convertFromCharPtr( const char* ptr )
    { return reinterpret_cast<const T*>(ptr); }
  // Indirect serialization
  static Ordinal fromCountToIndirectBytes(const Ordinal count, const T buffer[])
    { return fromCountToDirectBytes(count); }
  static void serialize(
    const Ordinal count, const T buffer[], const Ordinal bytes, char charBuffer[]
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(bytes!=fromCountToDirectBytes(count));
#else
      (void)count;
#endif
      const char *_buffer = convertToCharPtr(buffer);
      std::copy(_buffer,_buffer+bytes,charBuffer);
    }
  static Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                          const char charBuffer[])
    { return fromDirectBytesToCount(bytes); }
  static void deserialize(
    const Ordinal bytes, const char charBuffer[], const Ordinal count, T buffer[]
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(count!=fromDirectBytesToCount(bytes));
#endif
      char *_buffer = convertToCharPtr(buffer);
      std::copy(charBuffer,charBuffer+bytes,_buffer);
    }
};

// Whether 'char' is signed or unsigned depends on the implementation.
// However, on some systems (e.g., Clang 3.1 on Intel Mac), partially
// specializing for signed char and unsigned char, but not for char,
// does not work.  Thus, we include specializations for all three
// possibilities.
template<typename Ordinal>
class SerializationTraits<Ordinal,char>
  : public DirectSerializationTraits<Ordinal,char>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,signed char>
  : public DirectSerializationTraits<Ordinal,signed char>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned char>
  : public DirectSerializationTraits<Ordinal,unsigned char>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,short int>
  : public DirectSerializationTraits<Ordinal,short int>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned short int>
  : public DirectSerializationTraits<Ordinal,unsigned short int>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,int>
  : public DirectSerializationTraits<Ordinal,int>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned int>
  : public DirectSerializationTraits<Ordinal,unsigned int>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,long int>
  : public DirectSerializationTraits<Ordinal,long int>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned long int>
  : public DirectSerializationTraits<Ordinal,long unsigned int>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,float>
  : public DirectSerializationTraits<Ordinal,float>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,double>
  : public DirectSerializationTraits<Ordinal,double>
{};

#ifdef HAVE_TEUCHOS_LONG_DOUBLE
template<typename Ordinal>
class SerializationTraits<Ordinal,long double>
  : public DirectSerializationTraits<Ordinal,long double>
{};
#endif


// FIXME: How do we know that P1 and P2 are directly serializable?
template<typename Ordinal, typename P1, typename P2>
class SerializationTraits<Ordinal,std::pair<P1,P2> >
  : public DirectSerializationTraits<Ordinal,std::pair<P1,P2> >
{};

#ifdef HAVE_TEUCHOSCORE_QUADMATH
template<typename Ordinal>
class SerializationTraits<Ordinal,__float128>
  : public DirectSerializationTraits<Ordinal,__float128>
{};
#endif // HAVE_TEUCHOSCORE_QUADMATH

#ifdef HAVE_TEUCHOS_QD
template<typename Ordinal>
class SerializationTraits<Ordinal,dd_real>
  : public DirectSerializationTraits<Ordinal,dd_real>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,qd_real>
  : public DirectSerializationTraits<Ordinal,qd_real>
{};
#endif

#ifdef HAVE_TEUCHOS_COMPLEX

template<typename Ordinal>
class SerializationTraits<Ordinal,std::complex<float> >
  : public DirectSerializationTraits<Ordinal,std::complex<float> >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,std::complex<double> >
  : public DirectSerializationTraits<Ordinal,std::complex<double> >
{};

#endif // HAVE_TEUCHOS_COMPLEX

// Partial specialization for long long.
// On platforms with sizeof(ptrdiff_t) <= sizeof(long long),
// this should take care of the ptrdiff_t specialization as well,
// since we've covered all built-in signed integer types above
// with size <= sizeof(long long).
template<typename Ordinal>
class SerializationTraits<Ordinal, long long int>
  : public DirectSerializationTraits<Ordinal, long long int>
{};

// Partial specialization for unsigned long long.
// On platforms with sizeof(size_t) <= sizeof(unsigned long long),
// this should take care of the size_t specialization as well,
// since we've covered all built-in unsigned integer types above
// with size <= sizeof(unsigned long long).
template<typename Ordinal>
class SerializationTraits<Ordinal, unsigned long long int>
  : public DirectSerializationTraits<Ordinal, unsigned long long int>
{};

} // namespace Teuchos

#endif // TEUCHOS_SERIALIZATION_TRAITS_HPP
