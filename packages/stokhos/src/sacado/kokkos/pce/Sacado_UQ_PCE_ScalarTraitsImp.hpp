// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_UQ_PCE_SCALARTRAITSIMP_HPP
#define SACADO_UQ_PCE_SCALARTRAITSIMP_HPP

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_as.hpp"

#include <iterator>

namespace Sacado {

  namespace UQ {

    //! Implementation for Teuchos::ScalarTraits for all PCE types
    template <typename PCEType>
    struct PCEScalarTraitsImp {
      typedef typename PCEType::storage_type storage_type;
      typedef typename storage_type::value_type value_type;
      typedef typename storage_type::ordinal_type ordinal_type;
      typedef Teuchos::ScalarTraits<value_type> TVT;

      typedef typename TVT::magnitudeType value_mag_type;
      typedef typename TVT::halfPrecision value_half_type;
      typedef typename TVT::doublePrecision value_double_type;

      typedef typename Sacado::mpl::apply<storage_type,ordinal_type,value_mag_type>::type storage_mag_type;
      typedef typename Sacado::mpl::apply<storage_type,ordinal_type,value_half_type>::type storage_half_type;
      typedef typename Sacado::mpl::apply<storage_type,ordinal_type,value_double_type>::type storage_double_type;

      typedef value_mag_type magnitudeType;
      typedef typename Sacado::mpl::apply<PCEType, storage_half_type>::type halfPrecision;
      typedef typename Sacado::mpl::apply<PCEType, storage_double_type>::type doublePrecision;
      typedef typename TVT::coordinateType coordinateType;

      typedef value_type innerProductType;

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

      static magnitudeType magnitude(const PCEType& a) {
        return a.two_norm();
      }

      static innerProductType innerProduct(const PCEType& a, const PCEType& b) {
        return a.inner_product(b);
      }

      static PCEType zero()  { return PCEType(0.0); }

      static PCEType one()   { return PCEType(1.0); }

      static PCEType conjugate(const PCEType& x) {
        PCEType y = x;
        y.copyForWrite();
        y.val() = TVT::conjugate(x.val());
        return y;
      }

      static magnitudeType real(const PCEType& x) {
        magnitudeType m = magnitudeType(0.0);
        const ordinal_type sz = x.size();
        for (ordinal_type i=0; i<sz; ++i) {
          value_mag_type t = TVT::real(x.fastAccessCoeff(i));
          m +=t*t;
        }
        return std::sqrt(m);
      }

      static magnitudeType imag(const PCEType& x) {
        magnitudeType m = magnitudeType(0.0);
        const ordinal_type sz = x.size();
        for (ordinal_type i=0; i<sz; ++i) {
          value_mag_type t = TVT::imag(x.fastAccessCoeff(i));
          m +=t*t;
        }
        return std::sqrt(m);
      }


      static value_type nan() { return TVT::nan(); }

      static bool isnaninf(const PCEType& x) {
        for (int i=0; i<x.size(); i++)
          if (TVT::isnaninf(x.fastAccessCoeff(i)))
            return true;
        return false;
      }

      static void seedrandom(unsigned int s) { TVT::seedrandom(s); }

      static PCEType random() { return PCEType(TVT::random()); }

      static const char * name() { return "Sacado::UQ::PCE<>"; }

      static PCEType squareroot(const PCEType& x) { return std::sqrt(x); }

      static PCEType pow(const PCEType& x, const PCEType& y) {
        return std::pow(x,y);
      }

      static PCEType log(const PCEType& x) { return std::log(x); }

      static PCEType log10(const PCEType& x) { return std::log10(x); }

    }; // class PCEScalarTraitsImp

    //! Implementation for Teuchos::ValueTypeConversionTraits for all PCE types
    template <typename TypeTo, typename PCEType>
    struct PCEValueTypeConversionTraitsImp {
      typedef typename Sacado::ValueType<PCEType>::type ValueT;
      typedef Teuchos::ValueTypeConversionTraits<TypeTo,ValueT> VTCT;
      static TypeTo convert( const PCEType t ) {
        return VTCT::convert(t.val());
      }
      static TypeTo safeConvert( const PCEType t ) {
        return VTCT::safeConvert(t.val());
      }
    };

    //! Implementation of Teuchos::SerializationTraits for all PCE types
    template <typename Ordinal, typename PCEType>
    class PCESerializationTraitsImp {
      typedef typename Sacado::ValueType<PCEType>::type ValueT;
      typedef Teuchos::SerializationTraits<Ordinal,ValueT> vSerT;
      typedef Teuchos::SerializationTraits<Ordinal,int> iSerT;
      typedef Teuchos::SerializationTraits<Ordinal,Ordinal> oSerT;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization = false;

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      static Ordinal fromCountToIndirectBytes(const Ordinal count,
                                              const PCEType buffer[]) {
        Ordinal bytes = 0;
        for (Ordinal i=0; i<count; i++) {
          int sz = buffer[i].size();
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &sz);
          Ordinal b2 = vSerT::fromCountToIndirectBytes(sz, buffer[i].coeff());
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          bytes += b1+b2+b3;
        }
        return bytes;
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      static void serialize (const Ordinal count,
                             const PCEType buffer[],
                             const Ordinal bytes,
                             char charBuffer[]) {
        for (Ordinal i=0; i<count; i++) {
          // First serialize size
          int sz = buffer[i].size();
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &sz);
          iSerT::serialize(1, &sz, b1, charBuffer);
          charBuffer += b1;

          // Next serialize PCE coefficients
          Ordinal b2 = vSerT::fromCountToIndirectBytes(sz, buffer[i].coeff());
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          oSerT::serialize(1, &b2, b3, charBuffer);
          charBuffer += b3;
          vSerT::serialize(sz, buffer[i].coeff(), b2, charBuffer);
          charBuffer += b2;
        }
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      static Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                              const char charBuffer[]) {
        Ordinal count = 0;
        Ordinal bytes_used = 0;
        while (bytes_used < bytes) {

          // Bytes for size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          bytes_used += b1;
          charBuffer += b1;

          // Bytes for PCE coefficients
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
      static void deserialize (const Ordinal bytes,
                               const char charBuffer[],
                               const Ordinal count,
                               PCEType buffer[]) {
        for (Ordinal i=0; i<count; i++) {

          // Deserialize size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          const int *sz = iSerT::convertFromCharPtr(charBuffer);
          charBuffer += b1;

          // Make sure PCE object is ready to receive values
          // We assume it has already been initialized with the proper
          // cijk object
          if (buffer[i].size() != *sz)
            buffer[i].reset(buffer[i].cijk(), *sz);
          buffer[i].copyForWrite();

          // Deserialize PCE coefficients
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          charBuffer += b3;
          vSerT::deserialize(*b2, charBuffer, *sz, buffer[i].coeff());
          charBuffer += *b2;
        }

      }

      //@}

    };


    //! Serializer object for all PCE types
    template <typename Ordinal, typename PCEType, typename ValueSerializer>
    class PCESerializerImp {

    public:

      //! Typename of value serializer
      typedef ValueSerializer value_serializer_type;

      //! Typename of cijk
      typedef typename PCEType::cijk_type cijk_type;


    protected:
      typedef typename Sacado::ValueType<PCEType>::type ValueT;
      typedef Teuchos::SerializationTraits<Ordinal,int> iSerT;
      typedef Teuchos::SerializationTraits<Ordinal,Ordinal> oSerT;

      cijk_type cijk;
      Teuchos::RCP<const ValueSerializer> vs;
      int sz;

    public:

      /// \brief Whether the type T supports direct serialization.
      static const bool supportsDirectSerialization = false;

      PCESerializerImp(const cijk_type& cijk_,
                       const Teuchos::RCP<const ValueSerializer>& vs_) :
        cijk(cijk_), vs(vs_), sz(cijk.dimension()) {}

      //! Return specified serializer size
      cijk_type getSerializerCijk() const { return cijk; }

      //! Get nested value serializer
      Teuchos::RCP<const value_serializer_type> getValueSerializer() const {
        return vs; }

      //! @name Indirect serialization functions (always defined and supported)
      //@{

      /** \brief Return the number of bytes for <tt>count</tt> objects. */
      Ordinal fromCountToIndirectBytes(const Ordinal count,
                                       const PCEType buffer[]) const {
        Ordinal bytes = 0;
        PCEType *x = NULL;
        const PCEType *cx;
        for (Ordinal i=0; i<count; i++) {
          int my_sz = buffer[i].size();
          if (sz != my_sz) {
            if (x == NULL)
              x = new PCEType;
            *x = buffer[i];
            x->reset(cijk);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &sz);
          Ordinal b2 = vs->fromCountToIndirectBytes(sz, cx->coeff());
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          bytes += b1+b2+b3;
        }
        if (x != NULL)
          delete x;
        return bytes;
      }

      /** \brief Serialize to an indirect <tt>char[]</tt> buffer. */
      void serialize (const Ordinal count,
                      const PCEType buffer[],
                      const Ordinal bytes,
                      char charBuffer[]) const {
        PCEType *x = NULL;
        const PCEType *cx;
        for (Ordinal i=0; i<count; i++) {
          // First serialize size
          int my_sz = buffer[i].size();
          if (sz != my_sz) {
            if (x == NULL)
              x = new PCEType(cijk);
            *x = buffer[i];
            x->reset(cijk);
            cx = x;
          }
          else
            cx = &(buffer[i]);
          Ordinal b1 = iSerT::fromCountToIndirectBytes(1, &sz);
          iSerT::serialize(1, &sz, b1, charBuffer);
          charBuffer += b1;

          // Next serialize PCE coefficients
          Ordinal b2 = vs->fromCountToIndirectBytes(sz, cx->coeff());
          Ordinal b3 = oSerT::fromCountToIndirectBytes(1, &b2);
          oSerT::serialize(1, &b2, b3, charBuffer);
          charBuffer += b3;
          vs->serialize(sz, cx->coeff(), b2, charBuffer);
          charBuffer += b2;
        }
        if (x != NULL)
          delete x;
      }

      /** \brief Return the number of objects for <tt>bytes</tt> of storage. */
      Ordinal fromIndirectBytesToCount(const Ordinal bytes,
                                       const char charBuffer[]) const {
        Ordinal count = 0;
        Ordinal bytes_used = 0;
        while (bytes_used < bytes) {

          // Bytes for size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          bytes_used += b1;
          charBuffer += b1;

          // Bytes for PCE coefficients
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
      void deserialize (const Ordinal bytes,
                        const char charBuffer[],
                        const Ordinal count,
                        PCEType buffer[]) const {
        for (Ordinal i=0; i<count; i++) {

          // Deserialize size
          Ordinal b1 = iSerT::fromCountToDirectBytes(1);
          const int *my_sz = iSerT::convertFromCharPtr(charBuffer);
          charBuffer += b1;

          // Create empty PCE object of given size
          buffer[i].reset(cijk);

          // Deserialize PCE coefficients
          Ordinal b3 = oSerT::fromCountToDirectBytes(1);
          const Ordinal *b2 = oSerT::convertFromCharPtr(charBuffer);
          charBuffer += b3;
          vs->deserialize(*b2, charBuffer, *my_sz, buffer[i].coeff());
          charBuffer += *b2;
        }

      }

      //@}

    };

  } // namespace UQ

} // namespace Sacado

#endif // SACADO_FAD_SCALARTRAITSIMP_HPP
