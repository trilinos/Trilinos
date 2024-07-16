// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TAY_SERIALIZATIONTRAITSIMP_HPP
#define SACADO_TAY_SERIALIZATIONTRAITSIMP_HPP

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_SACADO_TEUCHOSCOMM

#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_RCP.hpp"

namespace Sacado {

  namespace Tay {

    //! Serialization implementation for all Taylor types
    template <typename Ordinal, typename TayType, typename Serializer>
    struct SerializationImp {

    private:

      //! How to serialize ints
      typedef Teuchos::SerializationTraits<Ordinal,unsigned int> iSerT;

      //! How to serialize ordinals
      typedef Teuchos::SerializationTraits<Ordinal,Ordinal> oSerT;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization = false;

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Serializer& vs,
                                              const Ordinal count,
                                              const TayType buffer[],
                                              const Ordinal sz = 0) {
        Ordinal bytes = 0;
        TayType *x = NULL;
        const TayType *cx;
        for (Ordinal i=0; i<count; i++) {
          unsigned int my_sz = buffer[i].degree()+1;
          unsigned int tot_sz = sz;
          if (sz == 0) tot_sz = my_sz;
          if (tot_sz != my_sz) {
            x = new TayType(buffer[i]);
            x->resize(tot_sz-1, true);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &tot_sz);
          Ordinal b2 = vs.fromCountToIndirectBytes(
            tot_sz, &(cx->fastAccessCoeff(0)));
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          bytes += b1+b2+b3;
          if (x != NULL) {
            delete x;
            x = NULL;
          }
        }
        return bytes;
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Serializer& vs,
                             const Ordinal count,
                             const TayType buffer[],
                             const Ordinal bytes,
                             char charBuffer[],
                             const Ordinal sz = 0) {
        TayType *x = NULL;
        const TayType *cx;
        for (Ordinal i=0; i<count; i++) {
          // First serialize degree
          unsigned int my_sz = buffer[i].degree()+1;
          unsigned int tot_sz = sz;
          if (sz == 0) tot_sz = my_sz;
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &tot_sz);
          iSerT::serialize(1, &tot_sz, b1, charBuffer);
          charBuffer += b1;

          // Next serialize taylor coefficients
          if (tot_sz != my_sz) {
            x = new TayType(buffer[i]);
            x->resize(tot_sz-1, true);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          Ordinal b2 = vs.fromCountToIndirectBytes(
            tot_sz, &(cx->fastAccessCoeff(0)));
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          oSerT::serialize(1, &b2, b3, charBuffer);
          charBuffer += b3;
          vs.serialize(tot_sz, &(cx->fastAccessCoeff(0)), b2, charBuffer);
          charBuffer += b2;
          if (x != NULL) {
            delete x;
            x = NULL;
          }
        }
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromIndirectBytesToCount(const Serializer& vs,
                                              const Ordinal bytes,
                                              const char charBuffer[],
                                              const Ordinal sz = 0) {
        Ordinal count = 0;
        Ordinal bytes_used = 0;
        while (bytes_used < bytes) {

          // Bytes for degree
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          bytes_used += b1;
          charBuffer += b1;

          // Bytes for taylor coefficients
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          bytes_used += b3;
          charBuffer += b3;
          bytes_used += *b2;
          charBuffer += *b2;

          ++count;
        }
        return count;
      }

      /** \brief Deserialize from an indirect <tt>char[]</tt> buffer. */
      static void deserialize (const Serializer& vs,
                               const Ordinal bytes,
                               const char charBuffer[],
                               const Ordinal count,
                               TayType buffer[],
                               const Ordinal sz = 0) {
        for (Ordinal i=0; i<count; i++) {

          // Deserialize degree
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          const unsigned int *my_sz = iSerT::convertFromCharPtr(charBuffer);
          charBuffer += b1;

          // Create empty Taylor object of given size
          unsigned int tot_sz = sz;
          if (sz == 0) tot_sz = *my_sz;
          buffer[i] = TayType(tot_sz-1, 0.0);

          // Deserialize taylor coefficients
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          charBuffer += b3;
          vs.deserialize(*b2, charBuffer, *my_sz,
                         &(buffer[i].fastAccessCoeff(0)));
          charBuffer += *b2;
        }

      }

      //@}

    };

    //! Implementation of Teuchos::SerializationTraits for all Taylor types
    template <typename Ordinal, typename TayType>
    struct SerializationTraitsImp {

    private:

      //! Value type of Taylor type
      typedef typename Sacado::ValueType<TayType>::type ValueT;

      //! Default serializer for values
      typedef Teuchos::DefaultSerializer<Ordinal,ValueT> DS;

      //! Default serializer type for values
      typedef typename DS::DefaultSerializerType ValueSerializer;

      //! Implementation
      typedef SerializationImp<Ordinal,TayType,ValueSerializer> Imp;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization =
        Imp::supportsDirectSerialization;

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                              const TayType buffer[]) {
        return Imp::fromCountToIndirectBytes(
          DS::getDefaultSerializer(), count, buffer);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Ordinal count,
                             const TayType buffer[],
                             const Ordinal bytes,
                             char charBuffer[]) {
        Imp::serialize(
          DS::getDefaultSerializer(), count, buffer, bytes, charBuffer);
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                              const char charBuffer[]) {
        return Imp::fromIndirectBytesToCount(
          DS::getDefaultSerializer(), bytes, charBuffer);
      }

      /** \brief Deserialize from an indirect <tt>char[]</tt> buffer. */
      static void deserialize (const Ordinal bytes,
                               const char charBuffer[],
                               const Ordinal count,
                               TayType buffer[]) {
        Imp::deserialize(
          DS::getDefaultSerializer(), bytes, charBuffer, count, buffer);
      }

      //@}

    };

    //! An indirect serialization object for all Taylor types
    template <typename Ordinal, typename TayType, typename ValueSerializer>
    class SerializerImp {

    private:

      //! Implementation
      typedef SerializationImp<Ordinal,TayType,ValueSerializer> Imp;

      //! Serializer for value types
      Teuchos::RCP<const ValueSerializer> vs;

      //! Specified number of derivative components;
      Ordinal sz;

    public:

      //! Typename of value serializer
      typedef ValueSerializer value_serializer_type;

      /// \brief Whether we support direct serialization.
      static const bool supportsDirectSerialization =
        Imp::supportsDirectSerialization;

      //! Constructor
      SerializerImp(const Teuchos::RCP<const ValueSerializer>& vs_,
                    Ordinal sz_ = 0) :
        vs(vs_), sz(sz_) {}

      //! Return specified serializer size
      Ordinal getSerializerSize() const { return sz; }

      //! Get nested value serializer
      Teuchos::RCP<const value_serializer_type> getValueSerializer() const {
        return vs; }

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      Ordinal fromCountToIndirectBytes(const Ordinal count,
                                       const TayType buffer[]) const {
        return Imp::fromCountToIndirectBytes(*vs, count, buffer, sz);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      void serialize (const Ordinal count,
                      const TayType buffer[],
                      const Ordinal bytes,
                      char charBuffer[]) const {
        Imp::serialize(*vs, count, buffer, bytes, charBuffer, sz);
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                       const char charBuffer[]) const {
        return Imp::fromIndirectBytesToCount(*vs, bytes, charBuffer, sz);
      }

      /** \brief Deserialize from an indirect <tt>char[]</tt> buffer. */
      void deserialize (const Ordinal bytes,
                        const char charBuffer[],
                        const Ordinal count,
                        TayType buffer[]) const {
        return Imp::deserialize(*vs, bytes, charBuffer, count, buffer, sz);
      }

      //@}

    };

  } // namespace Tay

} // namespace Sacado

#endif // HAVE_SACADO_TEUCHOSCOMM

#endif // SACADO_FAD_SERIALIZATIONTRAITSIMP_HPP
