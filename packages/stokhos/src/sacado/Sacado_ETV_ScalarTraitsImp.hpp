// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_ETV_SCALAR_TRAITS_IMP_HPP
#define SACADO_ETV_SCALAR_TRAITS_IMP_HPP

#ifdef HAVE_SACADO_TEUCHOS

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {
  namespace ETV {

    //! Specializtion of Teuchos::ScalarTraits for all Vector types
    template <typename VecType>
    struct ScalarTraitsImp {
      typedef VecType ScalarType;
      typedef typename Sacado::ValueType<ScalarType>::type ValueT;
    
      typedef typename Sacado::mpl::apply<ScalarType,typename Teuchos::ScalarTraits<ValueT>::magnitudeType>::type magnitudeType;
      typedef typename Sacado::mpl::apply<ScalarType,typename Teuchos::ScalarTraits<ValueT>::halfPrecision>::type halfPrecision;
      typedef typename Sacado::mpl::apply<ScalarType,typename Teuchos::ScalarTraits<ValueT>::doublePrecision>::type doublePrecision;
      
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
      static magnitudeType magnitude(const ScalarType& a) {
	return std::fabs(a); 
      }
      static ValueT zero()  { 
	return ValueT(0.0); 
      }
      static ValueT one()   { 
	return ValueT(1.0); 
      }
      
      // Conjugate is only defined for real derivative components
      static ScalarType conjugate(const ScalarType& x) {
	int sz = x.size();
	ScalarType y(sz, ValueT(0.0));
	for (int i=0; i<sz; i++)
	  y.fastAccessCoeff(i) = 
	    Teuchos::ScalarTraits<ValueT>::conjugate(x.fastAccessCoeff(i));
	return y;
      }   
      
      // Real part is only defined for real derivative components
      static ScalarType real(const ScalarType& x) { 
	int sz = x.size();
	ScalarType y(sz, ValueT(0.0));
	for (int i=0; i<sz; i++)
	  y.fastAccessCoeff(i) = 
	    Teuchos::ScalarTraits<ValueT>::real(x.fastAccessCoeff(i));
	return y;
      }
      
      // Imaginary part is only defined for real derivative components
      static ScalarType imag(const ScalarType& x) { 
	int sz = x.size();
	ScalarType y(sz, ValueT(0.0));
	for (int i=0; i<sz; i++)
	  y.fastAccessCoeff(i) = 
	    Teuchos::ScalarTraits<ValueT>::imag(x.fastAccessCoeff(i));
	return y;
      }
    
      static ValueT nan() {
	return Teuchos::ScalarTraits<ValueT>::nan(); 
      }
      static bool isnaninf(const ScalarType& x) { 
	for (int i=0; i<x.size(); i++)
	  if (Teuchos::ScalarTraits<ValueT>::isnaninf(x.fastAccessCoeff(i)))
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
	return Sacado::StringName<ScalarType>::eval(); 
      }
      static ScalarType squareroot(const ScalarType& x) {
	return std::sqrt(x); 
      }
      static ScalarType pow(const ScalarType& x, const ScalarType& y) { 
	return std::pow(x,y); 
      }

    }; // class ScalarTraitsImp


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
    template <typename Ordinal, typename VecType>
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
    struct StaticSerializationTraitsImp {
      typedef typename Sacado::ValueType<VecType>::type ValueT;
      typedef Teuchos::SerializationTraits<Ordinal,ValueT> vSerT;
      typedef Teuchos::DirectSerializationTraits<Ordinal,VecType> DSerT;
      typedef Sacado::ETV::SerializationTraitsImp<Ordinal,VecType> STI;

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

#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETV_VECTORTRAITS_HPP
