// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
