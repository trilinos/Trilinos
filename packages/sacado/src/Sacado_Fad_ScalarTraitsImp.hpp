// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_SCALARTRAITSIMP_HPP
#define SACADO_FAD_SCALARTRAITSIMP_HPP

#ifdef HAVE_SACADO_TEUCHOS

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Sacado_mpl_apply.hpp"

#include <iterator>

namespace Sacado {

  namespace Fad {

    //! Implementation for Teuchos::ScalarTraits for all Fad types
    template <typename FadType>
    struct ScalarTraitsImp {
      typedef typename Sacado::ValueType<FadType>::type ValueT;

      typedef typename mpl::apply<FadType,typename Teuchos::ScalarTraits<ValueT>::magnitudeType>::type magnitudeType;
      typedef typename mpl::apply<FadType,typename Teuchos::ScalarTraits<ValueT>::halfPrecision>::type halfPrecision;
      typedef typename mpl::apply<FadType,typename Teuchos::ScalarTraits<ValueT>::doublePrecision>::type doublePrecision;

      static const bool isComplex = Teuchos::ScalarTraits<ValueT>::isComplex;
      static const bool isOrdinal = Teuchos::ScalarTraits<ValueT>::isOrdinal;
      static const bool isComparable =
        Teuchos::ScalarTraits<ValueT>::isComparable;
      static const bool hasMachineParameters =
        Teuchos::ScalarTraits<ValueT>::hasMachineParameters;
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType eps() {
        return Teuchos::ScalarTraits<ValueT>::eps();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType sfmin() {
        return Teuchos::ScalarTraits<ValueT>::sfmin();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType base()  {
        return Teuchos::ScalarTraits<ValueT>::base();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType prec()  {
        return Teuchos::ScalarTraits<ValueT>::prec();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType t()     {
        return Teuchos::ScalarTraits<ValueT>::t();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType rnd()   {
        return Teuchos::ScalarTraits<ValueT>::rnd();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType emin()  {
        return Teuchos::ScalarTraits<ValueT>::emin();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType rmin()  {
        return Teuchos::ScalarTraits<ValueT>::rmin();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType emax()  {
        return Teuchos::ScalarTraits<ValueT>::emax();
      }
      static typename Teuchos::ScalarTraits<ValueT>::magnitudeType rmax()  {
        return Teuchos::ScalarTraits<ValueT>::rmax();
      }
      static magnitudeType magnitude(const FadType& a) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
          a, "Error, the input value to magnitude(...) a = " << a <<
          " can not be NaN!" );
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(a) == false, std::runtime_error,
                           "Complex magnitude is not a differentiable "
                           "function of complex inputs.");
#endif
        //return std::fabs(a);
        magnitudeType b(a.size(),
                        Teuchos::ScalarTraits<ValueT>::magnitude(a.val()));
        if (Teuchos::ScalarTraits<ValueT>::real(a.val()) >= 0)
          for (int i=0; i<a.size(); i++)
            b.fastAccessDx(i) =
              Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessDx(i));
        else
          for (int i=0; i<a.size(); i++)
            b.fastAccessDx(i) =
              -Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessDx(i));
        return b;
      }
      static ValueT zero()  {
        return ValueT(0.0);
      }
      static ValueT one()   {
        return ValueT(1.0);
      }

      // Conjugate is only defined for real derivative components
      static FadType conjugate(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(x) == false, std::runtime_error,
                           "Complex conjugate is not a differentiable "
                           "function of complex inputs.");
#endif
        FadType y = x;
        y.val() = Teuchos::ScalarTraits<ValueT>::conjugate(x.val());
        return y;
      }

      // Real part is only defined for real derivative components
      static FadType real(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(x) == false, std::runtime_error,
                           "Real component is not a differentiable "
                           "function of complex inputs.");
#endif
        FadType y = x;
        y.val() = Teuchos::ScalarTraits<ValueT>::real(x.val());
        return y;
      }

      // Imaginary part is only defined for real derivative components
      static FadType imag(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(is_fad_real(x) == false, std::runtime_error,
                           "Imaginary component is not a differentiable "
                           "function of complex inputs.");
#endif
        return FadType(Teuchos::ScalarTraits<ValueT>::imag(x.val()));
      }

      static ValueT nan() {
        return Teuchos::ScalarTraits<ValueT>::nan();
      }
      static bool isnaninf(const FadType& x) {
        if (Teuchos::ScalarTraits<ValueT>::isnaninf(x.val()))
          return true;
        for (int i=0; i<x.size(); i++)
          if (Teuchos::ScalarTraits<ValueT>::isnaninf(x.dx(i)))
            return true;
        return false;
      }
      static void seedrandom(unsigned int s) {
        Teuchos::ScalarTraits<ValueT>::seedrandom(s);
      }
      static ValueT random() {
        return Teuchos::ScalarTraits<ValueT>::random();
      }
      static std::string name() {
        return Sacado::StringName<FadType>::eval();
      }
      static FadType squareroot(const FadType& x) {
#ifdef TEUCHOS_DEBUG
        TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
          x, "Error, the input value to squareroot(...) a = " << x <<
          " can not be NaN!" );
#endif
        return std::sqrt(x);
      }
      static FadType pow(const FadType& x, const FadType& y) {
        return std::pow(x,y);
      }

      // Helper function to determine whether a complex value is real
      static bool is_complex_real(const ValueT& x) {
        return
          Teuchos::ScalarTraits<ValueT>::magnitude(x-Teuchos::ScalarTraits<ValueT>::real(x)) == 0;
      }

      // Helper function to determine whether a Fad type is real
      static bool is_fad_real(const FadType& x) {
        if (x.size() == 0)
          return true;
        if (Teuchos::ScalarTraits<ValueT>::isComplex) {
          if (!is_complex_real(x.val()))
            return false;
          for (int i=0; i<x.size(); i++)
            if (!is_complex_real(x.fastAccessDx(i)))
              return false;
        }
        return true;
      }

    }; // class ScalarTraitsImp

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

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOSCORE) && defined(HAVE_SACADO_TEUCHOSKOKKOSCOMM) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_TeuchosCommAdapters.hpp"
#include "Kokkos_View_Fad.hpp"

namespace Teuchos {

template<typename Ordinal,
         typename T, typename L, typename D, typename M>
void broadcast(const Comm<Ordinal>& comm,
               const int rootRank,
               const Ordinal count,
               const Kokkos::View<T,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad>& buffer)
{
  typedef Kokkos::View<T,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad> view_type;
  typename view_type::array_type array_buffer = buffer;
  Ordinal array_count = count * buffer.storage_size();
  broadcast( comm, rootRank, array_count, array_buffer );
}

template<typename Ordinal,
         typename T, typename L, typename D, typename M,
         typename Serializer>
void broadcast(const Comm<Ordinal>& comm,
               const Serializer& serializer,
               const int rootRank,
               const Ordinal count,
               const Kokkos::View<T,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad>& buffer)
{
  typedef Kokkos::View<T,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad> view_type;
  typename view_type::array_type array_buffer = buffer;
  Ordinal array_count = count * buffer.storage_size();
  broadcast( comm, *(serializer.getValueSerializer()), rootRank,
             array_count, array_buffer );
}

}

#endif

#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_FAD_SCALARTRAITSIMP_HPP
