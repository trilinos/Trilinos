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

#ifndef TEUCHOS_SERIALIZATION_TRAITS_HELPERS_HPP
#define TEUCHOS_SERIALIZATION_TRAITS_HELPERS_HPP

#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos {

/** \brief A class for instantiating a default serialization object.
 *
 * The serialization buffer classes below are generalized beyond using the
 * SerializationTraits to use a general serialization object.  This is to 
 * allow for more general types of serialization, e.g., when other data needs
 * to be used.
 */
template <typename Ordinal, typename T>
class DefaultSerializer {
public:

  //! Typename of default serializer
  typedef SerializationTraits<Ordinal,T> DefaultSerializerType;

  //! Return the default serializer
  static DefaultSerializerType& getDefaultSerializer() {
    if (defaultSerializer_ == Teuchos::null)
      defaultSerializer_ = Teuchos::rcp(new DefaultSerializerType);
    return *defaultSerializer_;
  }

  //! Return RCP to the default serializer
  static Teuchos::RCP<DefaultSerializerType> getDefaultSerializerRCP() {
    if (defaultSerializer_ == Teuchos::null)
      defaultSerializer_ = Teuchos::rcp(new DefaultSerializerType);
    return defaultSerializer_;
  }

private:

  //! RCP to the default serializer
  static Teuchos::RCP<DefaultSerializerType> defaultSerializer_;
};

template <typename Ordinal, typename T>
Teuchos::RCP<typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
DefaultSerializer<Ordinal,T>::defaultSerializer_ = Teuchos::null;

/** \brief Encapsulate how an array of non-const objects with value sematics
 * is serialized into a <tt>char[]</tt> array.
 *
 * Default version templated on bool indicating whether direct serialization
 * is supported.  The default version is empty with specializations below
 * for direct and indirect serialization.
 */
template <typename Ordinal, typename T, typename Serializer,
	  bool direct = Serializer::supportsDirectSerialization>
class ValueTypeSerializationBufferImp {};

/** \brief Encapsulate how an array of const objects with value sematics is
 * serialized into a <tt>const char[]</tt> array.
 *
 * Default version templated on bool indicating whether direct serialization
 * is supported.  The default version is empty with specializations below
 * for direct and indirect serialization.
 */
template <typename Ordinal, typename T, typename Serializer,
	  bool direct = Serializer::supportsDirectSerialization>
class ConstValueTypeSerializationBufferImp {};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Default version templated on bool indicating whether direct serialization
 * is supported.  The default version is empty with specializations below
 * for direct and indirect serialization.
 */
template <typename Ordinal, typename T, typename Serializer,
	  bool direct = Serializer::supportsDirectSerialization>
class ValueTypeDeserializationBufferImp {};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Default version templated on bool indicating whether direct serialization
 * is supported.  The default version is empty with specializations below
 * for direct and indirect serialization.
 */
template <typename Ordinal, typename T, typename Serializer,
	  bool direct = Serializer::supportsDirectSerialization>
class ConstValueTypeDeserializationBufferImp {};

/** \brief Encapsulate how an array of non-const objects with value sematics
 * is serialized into a <tt>char[]</tt> array.
 *
 * Specialization for direct serialization.
 */
template <typename Ordinal, typename T, typename Serializer>
class ValueTypeSerializationBufferImp<Ordinal,T,Serializer,true> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeSerializationBufferImp(
    const Ordinal count, T buffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ValueTypeSerializationBufferImp();
  /** \brief . */
  char* getCharBuffer() const;
  /** \brief . */
  Ordinal getBytes() const;
  /** \brief . */
  const ArrayView<char> getCharBufferView() const;
private:
  Ordinal    count_;
  T          *buffer_;
  Ordinal    bytes_;
  char       *charBuffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ValueTypeSerializationBufferImp();
  ValueTypeSerializationBufferImp(const ValueTypeSerializationBufferImp&);
  ValueTypeSerializationBufferImp& operator=(const ValueTypeSerializationBufferImp&);
};

/** \brief Encapsulate how an array of const objects with value sematics is
 * serialized into a <tt>const char[]</tt> array.
 *
 * Specialization for direct serialization.
 */
template <typename Ordinal, typename T, typename Serializer>
class ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,true> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeSerializationBufferImp(
    const Ordinal count, const T buffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ConstValueTypeSerializationBufferImp();
  /** \brief . */
  const char* getCharBuffer() const;
  /** \brief . */
  Ordinal getBytes() const;
  /** \brief . */
  const ArrayView<const char> getCharBufferView() const;
private:
  Ordinal    count_;
  const T    *buffer_;
  Ordinal    bytes_;
  const char *charBuffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ConstValueTypeSerializationBufferImp();
  ConstValueTypeSerializationBufferImp(const ConstValueTypeSerializationBufferImp&);
  ConstValueTypeSerializationBufferImp& operator=(const ConstValueTypeSerializationBufferImp&);
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Specialization for direct serialization.
 */
template <typename Ordinal, typename T, typename Serializer>
class ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeDeserializationBufferImp(
    const Ordinal bytes, char charBuffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ValueTypeDeserializationBufferImp();
  /** \brief . */
  T* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  Ordinal    bytes_;
  char       *charBuffer_;
  Ordinal    count_;
  T          *buffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ValueTypeDeserializationBufferImp();
  ValueTypeDeserializationBufferImp(const ValueTypeDeserializationBufferImp&);
  ValueTypeDeserializationBufferImp& operator=(const ValueTypeDeserializationBufferImp&);
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Specialization for direct serialization.
 */
template <typename Ordinal, typename T, typename Serializer>
class ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeDeserializationBufferImp(
    const Ordinal bytes, const char charBuffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ConstValueTypeDeserializationBufferImp();
  /** \brief . */
  const T* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  Ordinal    bytes_;
  const char *charBuffer_;
  Ordinal    count_;
  const T    *buffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ConstValueTypeDeserializationBufferImp();
  ConstValueTypeDeserializationBufferImp(const ConstValueTypeDeserializationBufferImp&);
  ConstValueTypeDeserializationBufferImp& operator=(const ConstValueTypeDeserializationBufferImp&);
};

/** \brief Encapsulate how an array of non-const objects with value sematics
 * is serialized into a <tt>char[]</tt> array.
 *
 * Specialization for indirect serialization
 */
template <typename Ordinal, typename T, typename Serializer>
class ValueTypeSerializationBufferImp<Ordinal,T,Serializer,false> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeSerializationBufferImp(
    const Ordinal count, T buffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ValueTypeSerializationBufferImp();
  /** \brief . */
  char* getCharBuffer() const;
  /** \brief . */
  Ordinal getBytes() const;
  /** \brief . */
  const ArrayView<char> getCharBufferView() const;
private:
  Ordinal    count_;
  T          *buffer_;
  Ordinal    bytes_;
  mutable Array<char> charBuffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ValueTypeSerializationBufferImp();
  ValueTypeSerializationBufferImp(const ValueTypeSerializationBufferImp&);
  ValueTypeSerializationBufferImp& operator=(const ValueTypeSerializationBufferImp&);
};

/** \brief Encapsulate how an array of const objects with value sematics is
 * serialized into a <tt>const char[]</tt> array.
 *
 * Specialization for indirect serialization
 */
template <typename Ordinal, typename T, typename Serializer>
class ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,false> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeSerializationBufferImp(
    const Ordinal count, const T buffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ConstValueTypeSerializationBufferImp();
  /** \brief . */
  const char* getCharBuffer() const;
  /** \brief . */
  Ordinal getBytes() const;
  /** \brief . */
  const ArrayView<const char> getCharBufferView() const;
private:
  Ordinal    count_;
  const T    *buffer_;
  Ordinal    bytes_;
  Array<char> charBuffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ConstValueTypeSerializationBufferImp();
  ConstValueTypeSerializationBufferImp(const ConstValueTypeSerializationBufferImp&);
  ConstValueTypeSerializationBufferImp& operator=(const ConstValueTypeSerializationBufferImp&);
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Specialization for indirect serialization
 */
template <typename Ordinal, typename T, typename Serializer>
class ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeDeserializationBufferImp(
    const Ordinal bytes, char charBuffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ValueTypeDeserializationBufferImp();
  /** \brief . */
  T* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  Ordinal    bytes_;
  char       *charBuffer_;
  Ordinal    count_;
  mutable Array<T>   buffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ValueTypeDeserializationBufferImp();
  ValueTypeDeserializationBufferImp(const ValueTypeDeserializationBufferImp&);
  ValueTypeDeserializationBufferImp& operator=(const ValueTypeDeserializationBufferImp&);
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Specialization for indirect serialization
 */
template <typename Ordinal, typename T, typename Serializer>
class ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false> {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeDeserializationBufferImp(
    const Ordinal bytes, const char charBuffer[], 
    const RCP<const Serializer>& serializer
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ConstValueTypeDeserializationBufferImp();
  /** \brief . */
  const T* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  Ordinal    bytes_;
  const char *charBuffer_;
  Ordinal    count_;
  Array<T>   buffer_;
  RCP<const Serializer> serializer_;
  // Not defined and not to be called
  ConstValueTypeDeserializationBufferImp();
  ConstValueTypeDeserializationBufferImp(const ConstValueTypeDeserializationBufferImp&);
  ConstValueTypeDeserializationBufferImp& operator=(const ConstValueTypeDeserializationBufferImp&);
};


/** \brief Encapsulate how an array of non-const objects with value sematics
 * is serialized into a <tt>char[]</tt> array.
 */
template <typename Ordinal, typename T,
	  typename Serializer = typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
class ValueTypeSerializationBuffer : 
    public ValueTypeSerializationBufferImp<Ordinal,T,Serializer> {
public:
  typedef ValueTypeSerializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeSerializationBuffer(
    const Ordinal count, T buffer[], 
    const RCP<const Serializer>& serializer
    ) : Base(count,buffer,serializer) {}
};

/** \brief Encapsulate how an array of const objects with value sematics is
 * serialized into a <tt>const char[]</tt> array.
 */
template <typename Ordinal, typename T,
	  typename Serializer = typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
class ConstValueTypeSerializationBuffer : 
    public ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer> {
public:
  typedef ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeSerializationBuffer(
    const Ordinal count, const T buffer[], 
    const RCP<const Serializer>& serializer
    ) : Base(count,buffer,serializer) {}
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 */
template <typename Ordinal, typename T,
	  typename Serializer = typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
class ValueTypeDeserializationBuffer : 
    public ValueTypeDeserializationBufferImp<Ordinal,T,Serializer> {
public:
  typedef ValueTypeDeserializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeDeserializationBuffer(
    const Ordinal bytes, char charBuffer[], 
    const RCP<const Serializer>& serializer
    ) : Base(bytes,charBuffer,serializer) {}
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 */
template <typename Ordinal, typename T,
	  typename Serializer = typename DefaultSerializer<Ordinal,T>::DefaultSerializerType>
class ConstValueTypeDeserializationBuffer : 
    public ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer> {
public:
  typedef ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeDeserializationBuffer(
    const Ordinal bytes, const char charBuffer[], 
    const RCP<const Serializer>& serializer
    ) : Base(bytes,charBuffer,serializer) {}
};

/** \brief Encapsulate how an array of non-const objects with value sematics
 * is serialized into a <tt>char[]</tt> array.
 *
 * Specialization for the default serializer object type with a default 
 * argument for the serializer object parameter.
 */
template <typename Ordinal, typename T>
class ValueTypeSerializationBuffer<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> : 
    public ValueTypeSerializationBufferImp<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> {
public:
  typedef DefaultSerializer<Ordinal,T> DS;  // work around for parsing bug in gcc 4.1-4.2
  typedef typename DS::DefaultSerializerType Serializer;
  typedef ValueTypeSerializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeSerializationBuffer(
    const Ordinal count, T buffer[], 
    const RCP<const Serializer>& serializer = DS::getDefaultSerializerRCP()
    ) : Base(count,buffer,serializer) {}
};

/** \brief Encapsulate how an array of const objects with value sematics is
 * serialized into a <tt>const char[]</tt> array.
 *
 * Specialization for the default serializer object type with a default 
 * argument for the serializer object parameter.
 */
template <typename Ordinal, typename T>
class ConstValueTypeSerializationBuffer<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> : 
    public ConstValueTypeSerializationBufferImp<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> {
public:
  typedef DefaultSerializer<Ordinal,T> DS;  // work around for parsing bug in gcc 4.1-4.2
  typedef typename DS::DefaultSerializerType Serializer;
  typedef ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeSerializationBuffer(
    const Ordinal count, const T buffer[], 
    const RCP<const Serializer>& serializer = DS::getDefaultSerializerRCP()
    ) : Base(count,buffer,serializer) {}
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Specialization for the default serializer object type with a default 
 * argument for the serializer object parameter.
 */
template <typename Ordinal, typename T>
class ValueTypeDeserializationBuffer<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> : 
    public ValueTypeDeserializationBufferImp<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> {
public:
  typedef DefaultSerializer<Ordinal,T> DS;  // work around for parsing bug in gcc 4.1-4.2
  typedef typename DS::DefaultSerializerType Serializer;
  typedef ValueTypeDeserializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeDeserializationBuffer(
    const Ordinal bytes, char charBuffer[], 
    const RCP<const Serializer>& serializer = DS::getDefaultSerializerRCP()
    ) : Base(bytes,charBuffer,serializer) {}
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 *
 * Specialization for the default serializer object type with a default 
 * argument for the serializer object parameter.
 */
template <typename Ordinal, typename T>
class ConstValueTypeDeserializationBuffer<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> : 
    public ConstValueTypeDeserializationBufferImp<Ordinal,T,typename DefaultSerializer<Ordinal,T>::DefaultSerializerType> {
public:
  typedef DefaultSerializer<Ordinal,T> DS;  // work around for parsing bug in gcc 4.1-4.2
  typedef typename DS::DefaultSerializerType Serializer;
  typedef ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer> Base;
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeDeserializationBuffer(
    const Ordinal bytes, const char charBuffer[], 
    const RCP<const Serializer>& serializer = DS::getDefaultSerializerRCP()
    ) : Base(bytes,charBuffer,serializer) {}
};

// /////////////////////////////////////
// Template implementations for direct serialization

//
// ValueTypeSerializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
ValueTypeSerializationBufferImp(
  const Ordinal count, T buffer[], const RCP<const Serializer>& serializer
  )
  :count_(count), buffer_(buffer), serializer_(serializer)
{
  bytes_ = serializer_->fromCountToDirectBytes(count_);
  charBuffer_ = serializer_->convertToCharPtr(buffer_);
}

template <typename Ordinal, typename T, typename Serializer>
ValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
~ValueTypeSerializationBufferImp()
{
  // There is nothing to do since the type uses direct serialization!
}

template <typename Ordinal, typename T, typename Serializer>
char* ValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
getCharBuffer() const
{
  return charBuffer_;
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
getBytes() const
{
  return bytes_;
}


template <typename Ordinal, typename T, typename Serializer>
const ArrayView<char>
ValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
getCharBufferView() const
{
  return arrayView(charBuffer_, bytes_);
}


//
// ConstValueTypeSerializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
ConstValueTypeSerializationBufferImp(
  const Ordinal count, const T buffer[], const RCP<const Serializer>& serializer
  )
  :count_(count), buffer_(buffer), serializer_(serializer)
{
  bytes_ = serializer_->fromCountToDirectBytes(count_);
  charBuffer_ = serializer_->convertToCharPtr(buffer_);
}

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
~ConstValueTypeSerializationBufferImp()
{
  // There is nothing to do since the type uses direct serialization!
}

template <typename Ordinal, typename T, typename Serializer>
const char* ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
getCharBuffer() const
{
  return charBuffer_;
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
getBytes() const
{
  return bytes_;
}

template <typename Ordinal, typename T, typename Serializer>
const ArrayView<const char>
ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,true>::
getCharBufferView() const
{
  return arrayView(charBuffer_, bytes_);
}

//
// ValueTypeDeserializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
ValueTypeDeserializationBufferImp(
  const Ordinal bytes, char charBuffer[], const RCP<const Serializer>& serializer
  )
  :bytes_(bytes), charBuffer_(charBuffer), serializer_(serializer)
{
  count_ = serializer_->fromDirectBytesToCount(bytes_);
  buffer_ = serializer_->convertFromCharPtr(charBuffer_);
}

template <typename Ordinal, typename T, typename Serializer>
ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
~ValueTypeDeserializationBufferImp()
{
  // There is nothing to do since the type uses direct serialization!
}

template <typename Ordinal, typename T, typename Serializer>
T* ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
getBuffer() const
{
  return buffer_;
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
getCount() const
{
  return count_;
}

//
// ConstValueTypeDeserializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
ConstValueTypeDeserializationBufferImp(
  const Ordinal bytes, const char charBuffer[], const RCP<const Serializer>& serializer
  )
  :bytes_(bytes), charBuffer_(charBuffer), serializer_(serializer)
{
  count_ = serializer_->fromDirectBytesToCount(bytes_);
  buffer_ = serializer_->convertFromCharPtr(charBuffer_);
}

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
~ConstValueTypeDeserializationBufferImp()
{
  // There is nothing to do since the type uses direct serialization!
}

template <typename Ordinal, typename T, typename Serializer>
const T* ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
getBuffer() const
{
  return buffer_;
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,true>::
getCount() const
{
  return count_;
}


// /////////////////////////////////////
// Template implementations for indirect specializations

//
// ValueTypeSerializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
ValueTypeSerializationBufferImp(
  const Ordinal count, T buffer[], const RCP<const Serializer>& serializer
  )
  :count_(count), buffer_(buffer), serializer_(serializer)
{
  bytes_ = serializer_->fromCountToIndirectBytes(count_, buffer_);
  charBuffer_.resize(bytes_);
  serializer_->serialize(count_, buffer_, bytes_, &charBuffer_[0]);
}

template <typename Ordinal, typename T, typename Serializer>
ValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
~ValueTypeSerializationBufferImp()
{
  serializer_->deserialize(bytes_, &charBuffer_[0], count_, buffer_);
}

template <typename Ordinal, typename T, typename Serializer>
char* ValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
getCharBuffer() const
{
  return &charBuffer_[0];
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
getBytes() const
{
  return bytes_;
}

template <typename Ordinal, typename T, typename Serializer>
const ArrayView<char>
ValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
getCharBufferView() const
{
  return charBuffer_.view(0, bytes_);
}


//
// ConstValueTypeSerializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
ConstValueTypeSerializationBufferImp(
  const Ordinal count, const T buffer[], const RCP<const Serializer>& serializer
  )
  :count_(count), buffer_(buffer), serializer_(serializer)
{
  bytes_ = serializer_->fromCountToIndirectBytes(count_, buffer_);
  charBuffer_.resize(bytes_);
  serializer_->serialize(count_, buffer_, bytes_, &charBuffer_[0]);
}

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
~ConstValueTypeSerializationBufferImp()
{
}

template <typename Ordinal, typename T, typename Serializer>
const char* ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
getCharBuffer() const
{
  return &charBuffer_[0];
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
getBytes() const
{
  return bytes_;
}

template <typename Ordinal, typename T, typename Serializer>
const ArrayView<const char>
ConstValueTypeSerializationBufferImp<Ordinal,T,Serializer,false>::
getCharBufferView() const
{
  return charBuffer_.view(0, bytes_);
}

//
// ValueTypeDeserializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
ValueTypeDeserializationBufferImp(
  const Ordinal bytes, char charBuffer[], const RCP<const Serializer>& serializer
  )
  :bytes_(bytes), charBuffer_(charBuffer), serializer_(serializer)
{
  count_ = serializer_->fromIndirectBytesToCount(bytes_, charBuffer_);
  buffer_.resize(count_);
  serializer_->deserialize(bytes_, charBuffer_, count_, &buffer_[0]);
}

template <typename Ordinal, typename T, typename Serializer>
ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
~ValueTypeDeserializationBufferImp()
{
  serializer_->serialize(count_, &buffer_[0], bytes_, charBuffer_);
}

template <typename Ordinal, typename T, typename Serializer>
T* ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
getBuffer() const
{
  return &buffer_[0];
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
getCount() const
{
  return count_;
}

//
// ConstValueTypeDeserializationBufferImp
//

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
ConstValueTypeDeserializationBufferImp(
  const Ordinal bytes, const char charBuffer[], const RCP<const Serializer>& serializer
  )
  :bytes_(bytes), charBuffer_(charBuffer), serializer_(serializer)
{
  count_ = serializer_->fromIndirectBytesToCount(bytes_, charBuffer_);
  buffer_.resize(count_);
  serializer_->deserialize(bytes_, charBuffer_, count_, &buffer_[0]);
}

template <typename Ordinal, typename T, typename Serializer>
ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
~ConstValueTypeDeserializationBufferImp()
{
}

template <typename Ordinal, typename T, typename Serializer>
const T* ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
getBuffer() const
{
  return &buffer_[0];
}

template <typename Ordinal, typename T, typename Serializer>
Ordinal ConstValueTypeDeserializationBufferImp<Ordinal,T,Serializer,false>::
getCount() const
{
  return count_;
}

} // namespace Teuchos

#endif // TEUCHOS_SERIALIZATION_TRAITS_HELPERS_HPP
