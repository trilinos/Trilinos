// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


// NOTE: Before adding specializations of ScalarTraits, make sure that they do not duplicate 
// specializations already present in PyTrilinos (see packages/PyTrilinos/src/Teuchos_Traits.i)


#ifndef TEUCHOS_SERIALIZATION_TRAITS_HPP
#define TEUCHOS_SERIALIZATION_TRAITS_HPP

#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_TEUCHOS_QD
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#endif

namespace Teuchos {


template<typename T>
struct UndefinedSerializationTraits {
  //! This function should not compile if there is an attempt to instantiate!
  static inline T notDefined() {return(T::this_type_is_missing_a_specialization());}
};


/** \brief Serialization traits class for types that use value semantics.
 *
 * There are one of two modes associated with serialization.
 * 
 * ToDo: Finish documenation!
 */
template <typename Ordinal, typename T>
class SerializationTraits {
public:
  
  //! @name Seialization type selection 
  //@{

  /** \brief Determines if the type supports direct serialization. */
  static const bool supportsDirectSerialization = false;

  //@}

  //! @name Direct serialization functions (not defined if supportsDirectSerialization==false) 
  //@{

  /** \brief Return the number of bytes for <tt>count</tt> objects. */
  static Ordinal fromCountToDirectBytes(const Ordinal count) { (void)count; UndefinedSerializationTraits<T>::notDefined(); return 0; }

  /** \brief Convert the pointer type to <tt>char*</tt>. */
  static char* convertToCharPtr( T* ptr ) { (void)ptr; UndefinedSerializationTraits<T>::notDefined(); return 0; }

  /** \brief Convert the pointer type to <tt>const char*</tt>. */
  static const char* convertToCharPtr( const T* ptr ) { (void)ptr; UndefinedSerializationTraits<T>::notDefined(); return 0; }

  /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
  static Ordinal fromDirectBytesToCount(const Ordinal bytes) { (void)bytes; UndefinedSerializationTraits<T>::notDefined(); return 0; }

  /** \brief Convert the pointer type from <tt>char*</tt>. */
  static T* convertFromCharPtr( char* ptr ) { (void)ptr; UndefinedSerializationTraits<T>::notDefined(); return 0; }

  /** \brief Convert the pointer type from <tt>char*</tt>. */
  static const T* convertFromCharPtr( const char* ptr ) { (void)ptr; UndefinedSerializationTraits<T>::notDefined(); return 0; }

  //@}

  //! @name Indirect serialization functions (always defined and supported) 
  //@{

  /** \brief Return the number of bytes for <tt>count</tt> objects. */
  static Ordinal fromCountToIndirectBytes(const Ordinal count) { (void)count; UndefinedSerializationTraits<T>::notDefined(); return 0; }

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
  static void serialize(
    const Ordinal count, const T buffer[], const Ordinal bytes, char charBuffer[]
    )
    { (void)count; (void)buffer; (void)bytes; (void)charBuffer; UndefinedSerializationTraits<T>::notDefined(); }

  /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
  static Ordinal fromIndirectBytesToCount(const Ordinal bytes) { (void)bytes; UndefinedSerializationTraits<T>::notDefined(); return 0; }

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
  static void deserialize(
    const Ordinal bytes, const char charBuffer[], const Ordinal count, T buffer[]
    )
    { (void)bytes; (void)charBuffer; (void)count; (void)buffer;
      UndefinedSerializationTraits<T>::notDefined(); }
  
  //@}

};


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
  static Ordinal fromCountToIndirectBytes(const Ordinal count)
    { return fromCountToDirectBytes(count); }
  static void serialize(
    const Ordinal count, const T buffer[], const Ordinal bytes, char charBuffer[]
    )
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPT(bytes!=fromCountToIndirectBytes(count));
#endif
      const char *_buffer = convertToCharPtr(buffer);
      std::copy(_buffer,_buffer+bytes,charBuffer);
    }
  static Ordinal fromIndirectBytesToCount(const Ordinal bytes) 
    { return fromDirectBytesToCount(bytes); }
  static void deserialize(
    const Ordinal bytes, const char charBuffer[], const Ordinal count, T buffer[]
    )
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPT(count!=fromIndirectBytesToCount(bytes));
#endif
      char *_buffer = convertToCharPtr(buffer);
      std::copy(charBuffer,charBuffer+bytes,_buffer);
    }
};

template<typename Ordinal>
class SerializationTraits<Ordinal,char>
  : public DirectSerializationTraits<Ordinal,char>
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,short int>
  : public DirectSerializationTraits<Ordinal,short int>
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
class SerializationTraits<Ordinal,long unsigned int>
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

template<typename Ordinal, typename P1, typename P2>
class SerializationTraits<Ordinal,std::pair<P1,P2> >
  : public DirectSerializationTraits<Ordinal,std::pair<P1,P2> >
{};

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

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

template<typename Ordinal>
class SerializationTraits<Ordinal, long long int>
  : public DirectSerializationTraits<Ordinal, long long int>
{};

#endif // HAVE_TEUCHOS_LONG_LONG_INT

} // namespace Teuchos

#endif // TEUCHOS_SERIALIZATION_TRAITS_HPP
