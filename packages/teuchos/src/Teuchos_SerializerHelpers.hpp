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

#ifndef TEUCHOS_SERIALIZER_HELPERS_HPP
#define TEUCHOS_SERIALIZER_HELPERS_HPP

#include "Teuchos_Serializer.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos {

/** \brief Encapsulate how an array of non-const objects with reference
 * sematics is serialized into a <tt>char[]</tt> array and deserialized again.
 */
template <typename Ordinal, typename T>
class ReferenceTypeSerializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ReferenceTypeSerializationBuffer(
    const Serializer<Ordinal,T> &serializer
    ,const Ordinal count, T*const buffer[]
    );
  /** \brief Deserialize from the interal <tt>char[]</tt> buffer back to the
   * original <tt>T*[]</tt> buffer.
   */
  ~ReferenceTypeSerializationBuffer();
  /** \brief . */
  char* getCharBuffer() const;
  /** \brief . */
  Ordinal getBytes() const;
private:
  const Serializer<Ordinal,T>  &serializer_;
  Ordinal                      count_;
  T*const                      *buffer_;
  Array<char>                  charBuffer_;
  // Not defined and not to be called
  ReferenceTypeSerializationBuffer();
  ReferenceTypeSerializationBuffer(const ReferenceTypeSerializationBuffer&);
  ReferenceTypeSerializationBuffer& operator=(const ReferenceTypeSerializationBuffer&);
};

/** \brief Encapsulate how an array of const objects with reference sematics
 * is serialized into a <tt>char[]</tt> array.
 */
template <typename Ordinal, typename T>
class ConstReferenceTypeSerializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>char[]</tt> buffer. */
  ConstReferenceTypeSerializationBuffer(
    const Serializer<Ordinal,T> &serializer
    ,const Ordinal count, const T*const buffer[]
    );
  /** \brief Free the internal <tt>char[]</tt> buffer (no data to be written
   * back).
   */
  ~ConstReferenceTypeSerializationBuffer();
  /** \brief . */
  const char* getCharBuffer() const;
  /** \brief . */
  Ordinal getBytes() const;
private:
  const Serializer<Ordinal,T>  &serializer_;
  Ordinal                      count_;
  const T*const                *buffer_;
  Ordinal                      bytes_;
  Array<char>                  charBuffer_;
  // Not defined and not to be called
  ConstReferenceTypeSerializationBuffer();
  ConstReferenceTypeSerializationBuffer(const ConstReferenceTypeSerializationBuffer&);
  ConstReferenceTypeSerializationBuffer& operator=(const ConstReferenceTypeSerializationBuffer&);
};

/** \brief Encapsulate how an array of non-const objects with reference
 * sematics is deserialized from a <tt>char[]</tt> array and then serialized
 * back into the <tt>char[]</tt> buffer again.
 */
template <typename Ordinal, typename T>
class ReferenceTypeDeserializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>T*[]</tt> buffer. */
  ReferenceTypeDeserializationBuffer(
    const Serializer<Ordinal,T> &serializer
    ,const Ordinal bytes, char charBuffer[]
    );
  /** \brief Reserialize back to the <tt>char[]</tt> buffer from the internal
   * <tt>T*[]</tt> buffer.
   */
  ~ReferenceTypeDeserializationBuffer();
  /** \brief . */
  T*const* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  typedef Array<RCP<T> >  buffer_ptr_t;
  typedef Array<T*>               buffer_t;
  const Serializer<Ordinal,T>  &serializer_;
  Ordinal                      bytes_;
  char                         *charBuffer_;
  buffer_ptr_t                 buffer_ptr_;
  buffer_t                     buffer_;
  // Not defined and not to be called
  ReferenceTypeDeserializationBuffer();
  ReferenceTypeDeserializationBuffer(const ReferenceTypeDeserializationBuffer&);
  ReferenceTypeDeserializationBuffer& operator=(const ReferenceTypeDeserializationBuffer&);
};

/** \brief Encapsulate how an array of onst objects with reference sematics is
 * deserialized from a <tt>char[]</tt> array with memory being automatically
 * freed at destruction time.
 */
template <typename Ordinal, typename T>
class ConstReferenceTypeDeserializationBuffer {
public:
  /** \brief Serialize to an internally stored <tt>T*[]</tt> buffer. */
  ConstReferenceTypeDeserializationBuffer(
    const Serializer<Ordinal,T> &serializer
    ,const Ordinal bytes, const char charBuffer[]
    );
  /** \brief Reserialize back to the <tt>char[]</tt> buffer from the internal
   * <tt>T*[]</tt> buffer.
   */
  ~ConstReferenceTypeDeserializationBuffer();
  /** \brief . */
  const T*const* getBuffer() const;
  /** \brief . */
  Ordinal getCount() const;
private:
  typedef Array<RCP<T> >  buffer_ptr_t;
  typedef Array<T*>               buffer_t;
  const Serializer<Ordinal,T>  &serializer_;
  Ordinal                      bytes_;
  const char                   *charBuffer_;
  buffer_ptr_t                 buffer_ptr_;
  buffer_t                     buffer_;
  // Not defined and not to be called
  ConstReferenceTypeDeserializationBuffer();
  ConstReferenceTypeDeserializationBuffer(const ConstReferenceTypeDeserializationBuffer&);
  ConstReferenceTypeDeserializationBuffer& operator=(const ConstReferenceTypeDeserializationBuffer&);
};

// /////////////////////////////////////
// Template implementations

//
// ReferenceTypeSerializationBuffer
//

template <typename Ordinal, typename T>
ReferenceTypeSerializationBuffer<Ordinal,T>::ReferenceTypeSerializationBuffer(
  const Serializer<Ordinal,T> &serializer
  ,const Ordinal count, T*const buffer[]
  )
  :serializer_(serializer), count_(count), buffer_(buffer)
{
  const Ordinal bytes = serializer_.getBufferSize(count_);
  charBuffer_.resize(bytes);
  serializer_.serialize(count_,buffer_,bytes,&charBuffer_[0]);
}

template <typename Ordinal, typename T>
ReferenceTypeSerializationBuffer<Ordinal,T>::~ReferenceTypeSerializationBuffer()
{
  serializer_.deserialize(charBuffer_.size(),&charBuffer_[0],count_,buffer_);
}

template <typename Ordinal, typename T>
char* ReferenceTypeSerializationBuffer<Ordinal,T>::getCharBuffer() const
{
  typedef ReferenceTypeSerializationBuffer<Ordinal,T>* this_ptr_t;
  return &(const_cast<this_ptr_t>(this)->charBuffer_)[0];
  // The above const_cast is a better alternative to declaring charBuffer_ to
  // be mutable, in my opinion.
}

template <typename Ordinal, typename T>
Ordinal ReferenceTypeSerializationBuffer<Ordinal,T>::getBytes() const
{
  return charBuffer_.size();
}

//
// ConstReferenceTypeSerializationBuffer
//

template <typename Ordinal, typename T>
ConstReferenceTypeSerializationBuffer<Ordinal,T>::ConstReferenceTypeSerializationBuffer(
  const Serializer<Ordinal,T> &serializer
  ,const Ordinal count, const T*const buffer[]
  )
  :serializer_(serializer), count_(count), buffer_(buffer)
{
  const Ordinal bytes = serializer_.getBufferSize(count_);
  charBuffer_.resize(bytes);
  serializer_.serialize(count_,buffer_,bytes,&charBuffer_[0]);
}

template <typename Ordinal, typename T>
ConstReferenceTypeSerializationBuffer<Ordinal,T>::~ConstReferenceTypeSerializationBuffer()
{
  // No need to copy back from the char[] buffer!
}

template <typename Ordinal, typename T>
const char* ConstReferenceTypeSerializationBuffer<Ordinal,T>::getCharBuffer() const
{
  return &charBuffer_[0];
}

template <typename Ordinal, typename T>
Ordinal ConstReferenceTypeSerializationBuffer<Ordinal,T>::getBytes() const
{
  return charBuffer_.size();
}

//
// ReferenceTypeDeserializationBuffer
//

template <typename Ordinal, typename T>
ReferenceTypeDeserializationBuffer<Ordinal,T>::ReferenceTypeDeserializationBuffer(
  const Serializer<Ordinal,T> &serializer
  ,const Ordinal bytes, char charBuffer[]
  )
  :serializer_(serializer),bytes_(bytes),charBuffer_(charBuffer)
{
  const Ordinal extent = serializer_.getBufferSize(1);
  const Ordinal count = bytes_ / extent;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( bytes_ % extent == 0 ) );
#endif
  buffer_ptr_.resize(count);
  buffer_.resize(count);
  for( int i = 0; i < count; ++i ) {
    buffer_ptr_[i] = serializer_.createObj();
    buffer_[i] = &*buffer_ptr_[i];
  }
  serializer_.deserialize(
    bytes_,charBuffer_,count,&buffer_[0]
    );
}

template <typename Ordinal, typename T>
ReferenceTypeDeserializationBuffer<Ordinal,T>::~ReferenceTypeDeserializationBuffer()
{
  serializer_.serialize(
    buffer_.size(),&buffer_[0],bytes_,charBuffer_
    );
}

template <typename Ordinal, typename T>
T*const* ReferenceTypeDeserializationBuffer<Ordinal,T>::getBuffer() const
{
  typedef ReferenceTypeDeserializationBuffer<Ordinal,T>* this_ptr_t;
  return &(const_cast<this_ptr_t>(this)->buffer_)[0];
  // The above const_cast is a better alternative to declaring buffer_ to be
  // mutable, in my opinion.
}

template <typename Ordinal, typename T>
Ordinal ReferenceTypeDeserializationBuffer<Ordinal,T>::getCount() const
{
  return buffer_.size();
}

//
// ConstReferenceTypeDeserializationBuffer
//

template <typename Ordinal, typename T>
ConstReferenceTypeDeserializationBuffer<Ordinal,T>::ConstReferenceTypeDeserializationBuffer(
  const Serializer<Ordinal,T> &serializer
  ,const Ordinal bytes, const char charBuffer[]
  )
  :serializer_(serializer),bytes_(bytes),charBuffer_(charBuffer)
{
  const Ordinal extent = serializer_.getBufferSize(1);
  const Ordinal count = bytes_ / extent;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( bytes_ % extent == 0 ) );
#endif
  buffer_ptr_.resize(count);
  buffer_.resize(count);
  for( int i = 0; i < count; ++i ) {
    buffer_ptr_[i] = serializer_.createObj();
    buffer_[i] = &*buffer_ptr_[i];
  }
  serializer_.deserialize(
    bytes_,charBuffer_,count,&buffer_[0]
    );
}

template <typename Ordinal, typename T>
ConstReferenceTypeDeserializationBuffer<Ordinal,T>::~ConstReferenceTypeDeserializationBuffer()
{
  // We don't need to serialized back into charBuffer_[] since it is constant!
}

template <typename Ordinal, typename T>
const T*const* ConstReferenceTypeDeserializationBuffer<Ordinal,T>::getBuffer() const
{
  return &buffer_[0];
}

template <typename Ordinal, typename T>
Ordinal ConstReferenceTypeDeserializationBuffer<Ordinal,T>::getCount() const
{
  return buffer_.size();
}

} // namespace Teuchos

#endif // TEUCHOS_SERIALIZER_HELPERS_HPP
