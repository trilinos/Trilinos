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

#ifndef TEUCHOS_SERIALIZER_HPP
#define TEUCHOS_SERIALIZER_HPP

#include "Teuchos_RCP.hpp"

namespace Teuchos {

/** \brief Strategy interface for the indirect serializing and deserializing
 * objects of a given type handled using reference semantics.
 *
 * This interface serializes and deserializes objects of type <tt>T</tt> to
 * and from independent <tt>char[]</tt> buffer arrays.  Direct serialization
 * (i.e. just using reinterpret casts) is not possible using this interface.
 */
template<typename Ordinal, typename T>
class Serializer {
public:

  /** \brief . */
  virtual ~Serializer() {}
  
  /** \brief Return an estimate for the maximum storage for <tt>count</tt>
   * objects to be serialized.
   */
  virtual Ordinal getBufferSize(const Ordinal count) const = 0;
  
  /** \brief Serialize an object to a <tt>char[]</tt> buffer.
   *
   * \param  count
   *           [in] Number of objects to be serialized.
   * \param  objs
   *           [in] Array (length <tt>count</tt>) for the objects to be serialized.
   * \param  bytes
   *           [in] Length of the buffer <tt>charBuffer[]</tt>
   * \param  charBuffer
   *           [out] Array (length <tt>bytes</tt>) that contains the serialized objects.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>count > 0</tt>
   * <li><tt>objs != NULL</tt>
   * <li><tt>charBuffer != NULL</tt>
   * <tt><tt>bytes == getBufferSize(count)</tt>
   * </ul>
   */
  virtual void serialize(
    const Ordinal          count
    ,const T * const       objs[]
    ,const Ordinal         bytes
    ,char                  charBuffer[]
    ) const = 0;

  // const ArrayView<const Ptr<const T> >& objs
  // std::ostream& oStream

  /** \brief Create an object of type <tt>T</tt> to be serialized into.
   */
  virtual RCP<T> createObj() const = 0;
  
  /** \brief Deserialize an object from a <tt>char[]</tt> buffer.
   *
   * \param  bytes
   *           [in] Length of the buffer <tt>charBuffer[]</tt>
   * \param  charBuffer
   *           [in] Array (length <tt>bytes</tt>) that contains the serialized objects.
   * \param  count
   *           [in] Number of objects to be deserialized.
   * \param  objs
   *           [out] Array (length <tt>count</tt>) for the deserialized objects.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>bytes > 0</tt>
   * <li><tt>objs != NULL</tt>
   * <li><tt>charBuffer != NULL</tt>
   * <tt><tt>bytes == getBufferSize(count)</tt>
   * </ul>
   */
  virtual void deserialize(
    const Ordinal         bytes
    ,const char           charBuffer[]
    ,const Ordinal        count
    ,T * const            objs[]
    ) const = 0;

  // const ArrayView<const Ptr<T> >& objs
  // std::istream& iStream

};

} // namespace Teuchos

#endif // TEUCHOS_SERIALIZER_HPP
