// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_MP_SCALAR_TRAITS_IMP_HPP
#define SACADO_MP_SCALAR_TRAITS_IMP_HPP

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {
  namespace MP {

    template <typename S, bool reduct_across_vector>
    struct ScalarTraitsImp {};

    // Implementation of Teuchos::ScalarTraits where reductions are taken
    // across the components of MP::Vector.  In this case magnitudeType is
    // a scalar
    template <typename S>
    struct ScalarTraitsImp<S,true> {
      typedef Sacado::MP::Vector<S> ScalarType;
      typedef typename S::value_type value_type;
      typedef typename S::ordinal_type ordinal_type;
      typedef Teuchos::ScalarTraits<value_type> TVT;

      typedef typename TVT::magnitudeType value_mag_type;
      typedef typename TVT::halfPrecision value_half_type;
      typedef typename TVT::doublePrecision value_double_type;

      typedef typename Sacado::mpl::apply<S,ordinal_type,value_mag_type>::type storage_mag_type;
      typedef typename Sacado::mpl::apply<S,ordinal_type,value_half_type>::type storage_half_type;
      typedef typename Sacado::mpl::apply<S,ordinal_type,value_double_type>::type storage_double_type;

      typedef value_mag_type magnitudeType;
      typedef Sacado::MP::Vector<storage_half_type> halfPrecision;
      typedef Sacado::MP::Vector<storage_double_type> doublePrecision;
      typedef typename Teuchos::ScalarTraits<value_type>::coordinateType coordinateType;

      static const bool isComplex = TVT::isComplex;
      static const bool isOrdinal = TVT::isOrdinal;
      static const bool isComparable = TVT::isComparable;
      static const bool hasMachineParameters = TVT::hasMachineParameters;

      static value_mag_type eps() { return TVT::eps(); }

      static value_mag_type sfmin() { return TVT::sfmin(); }

      static value_mag_type base()  { return TVT::base(); }

      static value_mag_type prec()  { return TVT::prec(); }

      static value_mag_type t()     { return TVT::t(); }

      static value_mag_type rnd()   { return TVT::rnd(); }

      static value_mag_type emin()  { return TVT::emin(); }

      static value_mag_type rmin()  { return TVT::rmin(); }

      static value_mag_type emax()  { return TVT::emax(); }

      static value_mag_type rmax()  { return TVT::rmax(); }

      static magnitudeType magnitude(const ScalarType& a) {
        magnitudeType m = magnitudeType(0.0);
        const ordinal_type sz = a.size();
        for (ordinal_type i=0; i<sz; ++i) {
          value_mag_type t = TVT::magnitude(a.fastAccessCoeff(i));
          m +=t*t;
        }
        return std::sqrt(m);
      }

      static ScalarType zero()  { return ScalarType(0.0); }

      static ScalarType one()   { return ScalarType(1.0); }


      static ScalarType conjugate(const ScalarType& x) {
        int sz = x.size();
        ScalarType y(sz, value_type(0.0));
        for (int i=0; i<sz; i++)
          y.fastAccessCoeff(i) = TVT::conjugate(x.fastAccessCoeff(i));
        return y;
      }


      static magnitudeType real(const ScalarType& x) {
        magnitudeType m = magnitudeType(0.0);
        const ordinal_type sz = x.size();
        for (ordinal_type i=0; i<sz; ++i) {
          value_mag_type t = TVT::real(x.fastAccessCoeff(i));
          m +=t*t;
        }
        return std::sqrt(m);
      }


      static magnitudeType imag(const ScalarType& x) {
        magnitudeType m = magnitudeType(0.0);
        const ordinal_type sz = x.size();
        for (ordinal_type i=0; i<sz; ++i) {
          value_mag_type t = TVT::imag(x.fastAccessCoeff(i));
          m +=t*t;
        }
        return std::sqrt(m);
      }

      static value_type nan() { return TVT::nan(); }

      static bool isnaninf(const ScalarType& x) {
        for (int i=0; i<x.size(); i++)
          if (TVT::isnaninf(x.fastAccessCoeff(i)))
            return true;
        return false;
      }

      static void seedrandom(unsigned int s) { TVT::seedrandom(s); }

      static ScalarType random() { return ScalarType(TVT::random()); }

      static const char * name() { return "Sacado::MP::Vector<>"; }

      static ScalarType squareroot(const ScalarType& x) { return std::sqrt(x); }

      static ScalarType pow(const ScalarType& x, const ScalarType& y) {
        return std::pow(x,y);
      }

      static ScalarType log(const ScalarType& x) { return std::log(x); }

      static ScalarType log10(const ScalarType& x) { return std::log10(x); }

    }; // class ScalarTraitsImp<S,true>

    // Implementation of Teuchos::ScalarTraits where reductions are not taken
    // across the components of MP::Vector.  In this case magnitudeType is
    // an MP::Vector
    template <typename S>
    struct ScalarTraitsImp<S,false> {
      typedef Sacado::MP::Vector<S> ScalarType;
      typedef typename S::value_type value_type;
      typedef typename S::ordinal_type ordinal_type;
      typedef Teuchos::ScalarTraits<value_type> TVT;

      typedef typename TVT::magnitudeType value_mag_type;
      typedef typename TVT::halfPrecision value_half_type;
      typedef typename TVT::doublePrecision value_double_type;

      typedef typename Sacado::mpl::apply<S,ordinal_type,value_mag_type>::type storage_mag_type;
      typedef typename Sacado::mpl::apply<S,ordinal_type,value_half_type>::type storage_half_type;
      typedef typename Sacado::mpl::apply<S,ordinal_type,value_double_type>::type storage_double_type;

      typedef Sacado::MP::Vector<storage_mag_type> magnitudeType;
      typedef Sacado::MP::Vector<storage_half_type> halfPrecision;
      typedef Sacado::MP::Vector<storage_double_type> doublePrecision;
      typedef typename Teuchos::ScalarTraits<value_type>::coordinateType coordinateType;

      static const bool isComplex = TVT::isComplex;
      static const bool isOrdinal = TVT::isOrdinal;
      static const bool isComparable = TVT::isComparable;
      static const bool hasMachineParameters = TVT::hasMachineParameters;

      static value_mag_type eps() { return TVT::eps(); }

      static value_mag_type sfmin() { return TVT::sfmin(); }

      static value_mag_type base()  { return TVT::base(); }

      static value_mag_type prec()  { return TVT::prec(); }

      static value_mag_type t()     { return TVT::t(); }

      static value_mag_type rnd()   { return TVT::rnd(); }

      static value_mag_type emin()  { return TVT::emin(); }

      static value_mag_type rmin()  { return TVT::rmin(); }

      static value_mag_type emax()  { return TVT::emax(); }

      static value_mag_type rmax()  { return TVT::rmax(); }

      static magnitudeType magnitude(const ScalarType& a) {
        return std::fabs(a);
      }

      static ScalarType zero()  { return ScalarType(0.0); }

      static ScalarType one()   { return ScalarType(1.0); }


      static ScalarType conjugate(const ScalarType& x) {
        int sz = x.size();
        ScalarType y(sz, value_type(0.0));
        for (int i=0; i<sz; i++)
          y.fastAccessCoeff(i) = TVT::conjugate(x.fastAccessCoeff(i));
        return y;
      }

      static ScalarType real(const ScalarType& x) {
        int sz = x.size();
        ScalarType y(sz, value_type(0.0));
        for (int i=0; i<sz; i++)
          y.fastAccessCoeff(i) = TVT::real(x.fastAccessCoeff(i));
        return y;
      }

      static ScalarType imag(const ScalarType& x) {
        int sz = x.size();
        ScalarType y(sz, value_type(0.0));
        for (int i=0; i<sz; i++)
          y.fastAccessCoeff(i) = TVT::imag(x.fastAccessCoeff(i));
        return y;
      }

      static value_type nan() { return TVT::nan(); }

      static bool isnaninf(const ScalarType& x) {
        for (int i=0; i<x.size(); i++)
          if (TVT::isnaninf(x.fastAccessCoeff(i)))
            return true;
        return false;
      }

      static void seedrandom(unsigned int s) { TVT::seedrandom(s); }

      static ScalarType random() { return ScalarType(TVT::random()); }

      static const char * name() { return "Sacado::MP::Vector<>"; }

      static ScalarType squareroot(const ScalarType& x) { return std::sqrt(x); }

      static ScalarType pow(const ScalarType& x, const ScalarType& y) {
        return std::pow(x,y);
      }

      static ScalarType log(const ScalarType& x) { return std::log(x); }

      static ScalarType log10(const ScalarType& x) { return std::log10(x); }

    }; // class ScalarTraitsImp<S,false>

    //! Serialization implementation for all Vector types
    template <typename Ordinal, typename VecType, typename Serializer>
    struct SerializationImp {

    private:

      //! Value type
      typedef typename Sacado::ValueType<VecType>::type ValueT;

      //! How to serialize ints
      typedef Teuchos::SerializationTraits<Ordinal,int> iSerT;

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
                                              const VecType buffer[],
                                              const Ordinal sz = 0) {
        Ordinal bytes = 0;
        VecType *x = NULL;
        const VecType *cx;
        for (Ordinal i=0; i<count; i++) {
          int my_sz = buffer[i].size();
          int tot_sz = sz;
          if (sz == 0) tot_sz = my_sz;
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &tot_sz);
          if (tot_sz != my_sz) {
            if (x == NULL)
              x = new VecType;
            *x = buffer[i];
            x->reset(tot_sz);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          Ordinal b2 = vs.fromCountToIndirectBytes(tot_sz, cx->coeff());
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          bytes += b1+b2+b3;
        }
        if (x != NULL)
          delete x;
        return bytes;
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Serializer& vs,
                             const Ordinal count,
                             const VecType buffer[],
                             const Ordinal bytes,
                             char charBuffer[],
                             const Ordinal sz = 0) {
        VecType *x = NULL;
        const VecType *cx;
        for (Ordinal i=0; i<count; i++) {
          // First serialize size
          int my_sz = buffer[i].size();
          int tot_sz = sz;
          if (sz == 0) tot_sz = my_sz;
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &tot_sz);
          iSerT::serialize(1, &tot_sz, b1, charBuffer);
          charBuffer += b1;

          // Next serialize vector coefficients
          if (tot_sz != my_sz) {
            if (x == NULL)
              x = new VecType;
            *x = buffer[i];
            x->reset(tot_sz);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          Ordinal b2 = vs.fromCountToIndirectBytes(tot_sz, cx->coeff());
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          oSerT::serialize(1, &b2, b3, charBuffer);
          charBuffer += b3;
          vs.serialize(tot_sz, cx->coeff(), b2, charBuffer);
          charBuffer += b2;
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

          // Bytes for vector coefficients
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
                               VecType buffer[],
                               const Ordinal sz = 0) {
        for (Ordinal i=0; i<count; i++) {

          // Deserialize size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          const int *my_sz = iSerT::convertFromCharPtr(charBuffer);
          charBuffer += b1;

          // Create empty Vector object of given size
          int tot_sz = sz;
          if (sz == 0) tot_sz = *my_sz;
          buffer[i] = VecType(tot_sz, ValueT(0.0));

          // Deserialize vector coefficients
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          charBuffer += b3;
          vs.deserialize(*b2, charBuffer, *my_sz, buffer[i].coeff());
          charBuffer += *b2;
        }

      }
      //@}

    };

    //! Implementation of Teuchos::SerializationTraits for all Vector types
    template <typename Ordinal, typename VecType, bool is_static = false>
    struct SerializationTraitsImp {

    private:

      //! Value type of Vec type
      typedef typename Sacado::ValueType<VecType>::type ValueT;

      //! Default serializer for values
      typedef Teuchos::DefaultSerializer<Ordinal,ValueT> DS;

      //! Default serializer type for values
      typedef typename DS::DefaultSerializerType ValueSerializer;

      //! Implementation
      typedef SerializationImp<Ordinal,VecType,ValueSerializer> Imp;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization =
        Imp::supportsDirectSerialization;

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                              const VecType buffer[]) {
        return Imp::fromCountToIndirectBytes(
          DS::getDefaultSerializer(), count, buffer);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Ordinal count,
                             const VecType buffer[],
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
                               VecType buffer[]) {
        Imp::deserialize(
          DS::getDefaultSerializer(), bytes, charBuffer, count, buffer);
      }

      //@}

    };

    //! Implementation of Teuchos::SerializationTraits for all static Vec types
    template <typename Ordinal, typename VecType>
    struct SerializationTraitsImp<Ordinal, VecType, true> {
      typedef typename Sacado::ValueType<VecType>::type ValueT;
      typedef Teuchos::SerializationTraits<Ordinal,ValueT> vSerT;
      typedef Teuchos::DirectSerializationTraits<Ordinal,VecType> DSerT;
      typedef Sacado::MP::SerializationTraitsImp<Ordinal,VecType> STI;

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
      static char* convertToCharPtr( VecType* ptr ) {
        return DSerT::convertToCharPtr(ptr);
      }

      /** \brief Convert the pointer type to <tt>const char*</tt>. */
      static const char* convertToCharPtr( const VecType* ptr ) {
        return DSerT::convertToCharPtr(ptr);
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromDirectBytesToCount(const Ordinal bytes) {
        return DSerT::fromDirectBytesToCount(bytes);
      }

      /** \brief Convert the pointer type from <tt>char*</tt>. */
      static VecType* convertFromCharPtr( char* ptr ) {
        return DSerT::convertFromCharPtr(ptr);
      }

      /** \brief Convert the pointer type from <tt>char*</tt>. */
      static const VecType* convertFromCharPtr( const char* ptr ) {
        return DSerT::convertFromCharPtr(ptr);
      }

      //@}

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                              const VecType buffer[]) {
        if (supportsDirectSerialization)
          return DSerT::fromCountToIndirectBytes(count, buffer);
        else
          return STI::fromCountToIndirectBytes(count, buffer);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Ordinal count,
                             const VecType buffer[],
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
                               VecType buffer[]) {
        if (supportsDirectSerialization)
          return DSerT::deserialize(bytes, charBuffer, count, buffer);
        else
          return STI::deserialize(bytes, charBuffer, count, buffer);
      }

      //@}

    };

    //! An indirect serialization object for all Vector types
    template <typename Ordinal, typename VecType, typename ValueSerializer>
    class SerializerImp {

    private:

      //! Implementation
      typedef SerializationImp<Ordinal,VecType,ValueSerializer> Imp;

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
                                       const VecType buffer[]) const {
        return Imp::fromCountToIndirectBytes(*vs, count, buffer, sz);
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      void serialize (const Ordinal count,
                      const VecType buffer[],
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
                        VecType buffer[]) const {
        return Imp::deserialize(*vs, bytes, charBuffer, count, buffer, sz);
      }

      //@}

    };

  }

}

#endif // SACADO_MP_SCALAR_TRAITS_IMP_HPP
