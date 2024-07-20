// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_SERIALIZATIONTRAITSIMP_HPP
#define SACADO_FAD_SERIALIZATIONTRAITSIMP_HPP

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_SACADO_TEUCHOSCOMM

#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_RCP.hpp"

namespace Sacado {

  namespace Fad {

    //! Serialization mplementation for all Fad types
    template <typename Ordinal, typename FadType, typename Serializer>
    struct SerializationImp {

    private:

      //! How to serialize ints
      typedef Teuchos::SerializationTraits<Ordinal,int> iSerT;

      //! How to serialize ordinals
      typedef Teuchos::SerializationTraits<Ordinal,Ordinal> oSerT;

      //! Value type
      typedef typename Sacado::ValueType<FadType>::type value_type;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization = false;

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Serializer& vs,
                                              const Ordinal count,
                                              const FadType buffer[],
                                              const Ordinal sz = 0) {
        Ordinal bytes = 0;
        FadType *x = NULL;
        const FadType *cx;
        for (Ordinal i=0; i<count; i++) {
          int my_sz = buffer[i].size();
          int tot_sz = sz;
          if (sz == 0) tot_sz = my_sz;
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &tot_sz);
          Ordinal b2 = vs.fromCountToIndirectBytes(1, &(buffer[i].val()));
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          Ordinal b4;
          if (tot_sz != my_sz) {
            if (x == NULL)
              x = new FadType(tot_sz, 0.0);
            *x = buffer[i];
            x->expand(tot_sz);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          b4 = vs.fromCountToIndirectBytes(tot_sz, cx->dx());
          Ordinal b5 = oSerT::fromCountToIndirectBytes(1, &b4);
          bytes += b1+b2+b3+b4+b5;
        }
        if (x != NULL)
          delete x;
        return bytes;
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Serializer& vs,
                             const Ordinal count,
                             const FadType buffer[],
                             const Ordinal bytes,
                             char charBuffer[],
                             const Ordinal sz = 0) {
        FadType *x = NULL;
        const FadType *cx;
        for (Ordinal i=0; i<count; i++) {
          // First serialize size
          int my_sz = buffer[i].size();
          int tot_sz = sz;
          if (sz == 0) tot_sz = my_sz;
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &tot_sz);
          iSerT::serialize(1, &tot_sz, b1, charBuffer);
          charBuffer += b1;

          // Next serialize value
          Ordinal b2 = vs.fromCountToIndirectBytes(1, &(buffer[i].val()));
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          oSerT::serialize(1, &b2, b3, charBuffer);
          charBuffer += b3;
          vs.serialize(1, &(buffer[i].val()), b2, charBuffer);
          charBuffer += b2;

          // Next serialize derivative array
          Ordinal b4;
          if (tot_sz != my_sz) {
            if (x == NULL)
              x = new FadType(tot_sz, 0.0);
            *x = buffer[i];
            x->expand(tot_sz);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          b4 = vs.fromCountToIndirectBytes(tot_sz, cx->dx());
          Ordinal b5 = oSerT::fromCountToIndirectBytes(1, &b4);
          oSerT::serialize(1, &b4, b5, charBuffer);
          charBuffer += b5;
          vs.serialize(tot_sz, cx->dx(), b4, charBuffer);
          charBuffer += b4;
        }
        if (x != NULL)
          delete x;
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromIndirectBytesToCount(const Serializer& vs,
                                              const Ordinal bytes,
                                              const char charBuffer[],
                                              const Ordinal sz = 0) {
        Ordinal count = 0;
        Ordinal bytes_used = 0;
        while (bytes_used < bytes) {

          // Bytes for size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          bytes_used += b1;
          charBuffer += b1;

          // Bytes for value
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          bytes_used += b3;
          charBuffer += b3;
          bytes_used += *b2;
          charBuffer += *b2;

          // Bytes for derivative array
          Ordinal b5 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b4 = oSerT::convertFromCharPtr(charBuffer);
          bytes_used += b5;
          charBuffer += b5;
          bytes_used += *b4;
          charBuffer += *b4;

          ++count;
        }
        return count;
      }

      /** \brief Deserialize from an indirect <tt>char[]</tt> buffer. */
      static void deserialize (const Serializer& vs,
                               const Ordinal bytes,
                               const char charBuffer[],
                               const Ordinal count,
                               FadType buffer[],
                               const Ordinal sz = 0) {
        for (Ordinal i=0; i<count; i++) {

          // Deserialize size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          const int *my_sz = iSerT::convertFromCharPtr(charBuffer);
          charBuffer += b1;

          // Create empty Fad object of given size
          int tot_sz = sz;
          if (sz == 0) tot_sz = *my_sz;
          buffer[i] = FadType(tot_sz, 0.0);

          // Deserialize value
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          charBuffer += b3;
          vs.deserialize(*b2, charBuffer, 1, &(buffer[i].val()));
          charBuffer += *b2;

          // Deserialize derivative array
          Ordinal b5 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b4 = oSerT::convertFromCharPtr(charBuffer);
          charBuffer += b5;
          vs.deserialize(*b4, charBuffer, *my_sz,
                             &(buffer[i].fastAccessDx(0)));
          charBuffer += *b4;
        }

      }

      //@}

    };

    //! Implementation of Teuchos::SerializationTraits for all Fad types
    template <typename Ordinal, typename FadType>
    struct SerializationTraitsImp {

    private:

      //! Value type of Fad type
      typedef typename Sacado::ValueType<FadType>::type ValueT;

      //! Default serializer for values
      typedef Teuchos::DefaultSerializer<Ordinal,ValueT> DS;

      //! Default serializer type for values
      typedef typename DS::DefaultSerializerType ValueSerializer;

      //! Implementation
      typedef SerializationImp<Ordinal,FadType,ValueSerializer> Imp;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization =
        Imp::supportsDirectSerialization;

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                              const FadType buffer[]) {
        return Imp::fromCountToIndirectBytes(
          DS::getDefaultSerializer(), count, buffer);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Ordinal count,
                             const FadType buffer[],
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
                               FadType buffer[]) {
        Imp::deserialize(
          DS::getDefaultSerializer(), bytes, charBuffer, count, buffer);
      }

      //@}

    };

    //! Implementation of Teuchos::SerializationTraits for all static Fad types
    template <typename Ordinal, typename FadType>
    struct StaticSerializationTraitsImp {
      typedef typename Sacado::ValueType<FadType>::type ValueT;
      typedef Teuchos::SerializationTraits<Ordinal,ValueT> vSerT;
      typedef Teuchos::DirectSerializationTraits<Ordinal,FadType> DSerT;
      typedef Sacado::Fad::SerializationTraitsImp<Ordinal,FadType> STI;

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization =
        vSerT::supportsDirectSerialization;

      //! @name Direct serialization functions (not defined if supportsDirectSerialization==false)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToDirectBytes(const Ordinal count) {
        return DSerT::fromCountToDirectBytes(count);
      }

      /** \brief Convert the pointer type to <tt>char*</tt>. */
      static char* convertToCharPtr( FadType* ptr ) {
        return DSerT::convertToCharPtr(ptr);
      }

      /** \brief Convert the pointer type to <tt>const char*</tt>. */
      static const char* convertToCharPtr( const FadType* ptr ) {
        return DSerT::convertToCharPtr(ptr);
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromDirectBytesToCount(const Ordinal bytes) {
        return DSerT::fromDirectBytesToCount(bytes);
      }

      /** \brief Convert the pointer type from <tt>char*</tt>. */
      static FadType* convertFromCharPtr( char* ptr ) {
        return DSerT::convertFromCharPtr(ptr);
      }

      /** \brief Convert the pointer type from <tt>char*</tt>. */
      static const FadType* convertFromCharPtr( const char* ptr ) {
        return DSerT::convertFromCharPtr(ptr);
      }

      //@}

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                              const FadType buffer[]) {
        if (supportsDirectSerialization)
          return DSerT::fromCountToIndirectBytes(count, buffer);
        else
          return STI::fromCountToIndirectBytes(count, buffer);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Ordinal count,
                             const FadType buffer[],
                             const Ordinal bytes,
                             char charBuffer[]) {
        if (supportsDirectSerialization)
          return DSerT::serialize(count, buffer, bytes, charBuffer);
        else
          return STI::serialize(count, buffer, bytes, charBuffer);
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                              const char charBuffer[]) {
        if (supportsDirectSerialization)
          return DSerT::fromIndirectBytesToCount(bytes, charBuffer);
        else
          return STI::fromIndirectBytesToCount(bytes, charBuffer);
      }

      /** \brief Deserialize from an indirect <tt>char[]</tt> buffer. */
      static void deserialize (const Ordinal bytes,
                               const char charBuffer[],
                               const Ordinal count,
                               FadType buffer[]) {
        if (supportsDirectSerialization)
          return DSerT::deserialize(bytes, charBuffer, count, buffer);
        else
          return STI::deserialize(bytes, charBuffer, count, buffer);
      }

      //@}

    };

    //! An indirect serialization object for all Fad types
    template <typename Ordinal, typename FadType, typename ValueSerializer>
    class SerializerImp {

    private:

      //! Implementation
      typedef SerializationImp<Ordinal,FadType,ValueSerializer> Imp;

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
                                       const FadType buffer[]) const {
        return Imp::fromCountToIndirectBytes(*vs, count, buffer, sz);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      void serialize (const Ordinal count,
                      const FadType buffer[],
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
                        FadType buffer[]) const {
        return Imp::deserialize(*vs, bytes, charBuffer, count, buffer, sz);
      }

      //@}

    };

  } // namespace Fad

} // namespace Sacado

#endif // HAVE_SACADO_TEUCHOSCOMM

#endif // SACADO_FAD_SERIALIZATIONTRAITSIMP_HPP
