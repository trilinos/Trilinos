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

namespace Teuchos {

/** \brief Encapsulate how an array of non-const objects with value sematics
 * is serialized into a <tt>char[]</tt> array.
 */
template <typename Ordinal, typename T>
class ValueTypeSerializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeSerializationBuffer(
    const Ordinal count, T buffer[] 
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ValueTypeSerializationBuffer();
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
  // Not defined and not to be called
  ValueTypeSerializationBuffer();
  ValueTypeSerializationBuffer(const ValueTypeSerializationBuffer&);
  ValueTypeSerializationBuffer& operator=(const ValueTypeSerializationBuffer&);
};

/** \brief Encapsulate how an array of const objects with value sematics is
 * serialized into a <tt>const char[]</tt> array.
 */
template <typename Ordinal, typename T>
class ConstValueTypeSerializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeSerializationBuffer(
    const Ordinal count, const T buffer[]
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ConstValueTypeSerializationBuffer();
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
  // Not defined and not to be called
  ConstValueTypeSerializationBuffer();
  ConstValueTypeSerializationBuffer(const ConstValueTypeSerializationBuffer&);
  ConstValueTypeSerializationBuffer& operator=(const ConstValueTypeSerializationBuffer&);
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 */
template <typename Ordinal, typename T>
class ValueTypeDeserializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ValueTypeDeserializationBuffer(
    const Ordinal bytes, char charBuffer[] 
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ValueTypeDeserializationBuffer();
  /** \brief . */
  T* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  Ordinal    bytes_;
  char       *charBuffer_;
  Ordinal    count_;
  T          *buffer_;
  // Not defined and not to be called
  ValueTypeDeserializationBuffer();
  ValueTypeDeserializationBuffer(const ValueTypeDeserializationBuffer&);
  ValueTypeDeserializationBuffer& operator=(const ValueTypeDeserializationBuffer&);
};

/** \brief Encapsulate how an array of non-const serialized objects with value
 * sematics stored in a <tt>char[]</tt> array is deserialized to a
 * <tt>T[]</tt> array and then serialized back again.
 */
template <typename Ordinal, typename T>
class ConstValueTypeDeserializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstValueTypeDeserializationBuffer(
    const Ordinal bytes, const char charBuffer[]
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T[]</tt> buffer.
   */
  ~ConstValueTypeDeserializationBuffer();
  /** \brief . */
  const T* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  Ordinal    bytes_;
  const char *charBuffer_;
  Ordinal    count_;
  const T    *buffer_;
  // Not defined and not to be called
  ConstValueTypeDeserializationBuffer();
  ConstValueTypeDeserializationBuffer(const ConstValueTypeDeserializationBuffer&);
  ConstValueTypeDeserializationBuffer& operator=(const ConstValueTypeDeserializationBuffer&);
};

// /////////////////////////////////////
// Template implementations

//
// ValueTypeSerializationBuffer
//
// ToDo: Update this implementation to handle objects with indirect
// serialization when needed!
//

template <typename Ordinal, typename T>
ValueTypeSerializationBuffer<Ordinal,T>::ValueTypeSerializationBuffer(
  const Ordinal count, T buffer[]
  )
  :count_(count), buffer_(buffer)
{
  typedef SerializationTraits<Ordinal,T> SerT;
  bytes_ = SerT::fromCountToDirectBytes(count_);
  charBuffer_ = SerT::convertToCharPtr(buffer_);
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
ValueTypeSerializationBuffer<Ordinal,T>::~ValueTypeSerializationBuffer()
{
  // There is nothing to do since the type uses direct serialization!
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
char* ValueTypeSerializationBuffer<Ordinal,T>::getCharBuffer() const
{
  return charBuffer_;
}

template <typename Ordinal, typename T>
Ordinal ValueTypeSerializationBuffer<Ordinal,T>::getBytes() const
{
  return bytes_;
}


template <typename Ordinal, typename T>
const ArrayView<char>
ValueTypeSerializationBuffer<Ordinal,T>::getCharBufferView() const
{
  return arrayView(charBuffer_, bytes_);
}


//
// ConstValueTypeSerializationBuffer
//
// ToDo: Update this implementation to handle objects with indirect
// serialization when needed!
//

template <typename Ordinal, typename T>
ConstValueTypeSerializationBuffer<Ordinal,T>::ConstValueTypeSerializationBuffer(
  const Ordinal count, const T buffer[]
  )
  :count_(count), buffer_(buffer)
{
  typedef SerializationTraits<Ordinal,T> SerT;
  bytes_ = SerT::fromCountToDirectBytes(count_);
  charBuffer_ = SerT::convertToCharPtr(buffer_);
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
ConstValueTypeSerializationBuffer<Ordinal,T>::~ConstValueTypeSerializationBuffer()
{
  // There is nothing to do since the type uses direct serialization!
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
const char* ConstValueTypeSerializationBuffer<Ordinal,T>::getCharBuffer() const
{
  return charBuffer_;
}

template <typename Ordinal, typename T>
Ordinal ConstValueTypeSerializationBuffer<Ordinal,T>::getBytes() const
{
  return bytes_;
}

template <typename Ordinal, typename T>
const ArrayView<const char>
ConstValueTypeSerializationBuffer<Ordinal,T>::getCharBufferView() const
{
  return arrayView(charBuffer_, bytes_);
}

//
// ValueTypeDeserializationBuffer
//
// ToDo: Update this implementation to handle objects with indirect
// serialization when needed!
//

template <typename Ordinal, typename T>
ValueTypeDeserializationBuffer<Ordinal,T>::ValueTypeDeserializationBuffer(
  const Ordinal bytes, char charBuffer[]
  )
  :bytes_(bytes), charBuffer_(charBuffer)
{
  typedef SerializationTraits<Ordinal,T> SerT;
  count_ = SerT::fromDirectBytesToCount(bytes_);
  buffer_ = SerT::convertFromCharPtr(charBuffer_);
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
ValueTypeDeserializationBuffer<Ordinal,T>::~ValueTypeDeserializationBuffer()
{
  // There is nothing to do since the type uses direct serialization!
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
T* ValueTypeDeserializationBuffer<Ordinal,T>::getBuffer() const
{
  return buffer_;
}

template <typename Ordinal, typename T>
Ordinal ValueTypeDeserializationBuffer<Ordinal,T>::getCount() const
{
  return count_;
}

//
// ConstValueTypeDeserializationBuffer
//
// ToDo: Update this implementation to handle objects with indirect
// serialization when needed!
//

template <typename Ordinal, typename T>
ConstValueTypeDeserializationBuffer<Ordinal,T>::ConstValueTypeDeserializationBuffer(
  const Ordinal bytes, const char charBuffer[]
  )
  :bytes_(bytes), charBuffer_(charBuffer)
{
  typedef SerializationTraits<Ordinal,T> SerT;
  count_ = SerT::fromDirectBytesToCount(bytes_);
  buffer_ = SerT::convertFromCharPtr(charBuffer_);
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
ConstValueTypeDeserializationBuffer<Ordinal,T>::~ConstValueTypeDeserializationBuffer()
{
  // There is nothing to do since the type uses direct serialization!
  // ToDo: Handle indirect serailization!
}

template <typename Ordinal, typename T>
const T* ConstValueTypeDeserializationBuffer<Ordinal,T>::getBuffer() const
{
  return buffer_;
}

template <typename Ordinal, typename T>
Ordinal ConstValueTypeDeserializationBuffer<Ordinal,T>::getCount() const
{
  return count_;
}

} // namespace Teuchos

#endif // TEUCHOS_SERIALIZATION_TRAITS_HELPERS_HPP
