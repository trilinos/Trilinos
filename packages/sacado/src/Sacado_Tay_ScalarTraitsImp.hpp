// $Id$ 
// $Source$ 
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

#ifndef SACADO_TAY_SCALARTRAITSIMP_HPP
#define SACADO_TAY_SCALARTRAITSIMP_HPP

#ifdef HAVE_SACADO_TEUCHOS

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_apply.hpp"

#include <iterator>

namespace Sacado {

  namespace Tay {

    //! Implementation for Teuchos::ScalarTraits for all Taylor types
    template <typename TayType>
    struct ScalarTraitsImp {
      typedef typename Sacado::ValueType<TayType>::type ValueT;

      typedef typename mpl::apply<TayType,typename Teuchos::ScalarTraits<ValueT>::magnitudeType>::type magnitudeType;
      typedef typename mpl::apply<TayType,typename Teuchos::ScalarTraits<ValueT>::halfPrecision>::type halfPrecision;
      typedef typename mpl::apply<TayType,typename Teuchos::ScalarTraits<ValueT>::doublePrecision>::type doublePrecision;

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
      static magnitudeType magnitude(const TayType& a) {
#ifdef TEUCHOS_DEBUG
	TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
	  a, "Error, the input value to magnitude(...) a = " << a << 
	  " can not be NaN!" );
	TEUCHOS_TEST_FOR_EXCEPTION(is_tay_real(a) == false, std::runtime_error,
			   "Complex magnitude is not a differentiable "
			   "function of complex inputs.");
#endif
	//return std::fabs(a); 
	magnitudeType b(a.degree(), 
			Teuchos::ScalarTraits<ValueT>::magnitude(a.val()));
	if (Teuchos::ScalarTraits<ValueT>::real(a.val()) >= 0)
	  for (int i=1; i<=a.degree(); i++)
	    b.fastAccessCoeff(i) = 
	      Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessCoeff(i));
	else
	  for (int i=1; i<=a.degree(); i++)
	    b.fastAccessCoeff(i) = 
	      -Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessCoeff(i));
	return b;
      }
      static ValueT zero()  { 
	return ValueT(0.0); 
      }
      static ValueT one()   { 
	return ValueT(1.0); 
      }
      
      // Conjugate is only defined for real derivative components
      static TayType conjugate(const TayType& x) {
#ifdef TEUCHOS_DEBUG
	TEUCHOS_TEST_FOR_EXCEPTION(is_tay_real(x) == false, std::runtime_error,
			   "Complex conjugate is not a differentiable "
			   "function of complex inputs.");
#endif
	TayType y = x;
	y.copyForWrite();
	y.val() = Teuchos::ScalarTraits<ValueT>::conjugate(x.val());
	return y;
      }   

      // Real part is only defined for real derivative components
      static TayType real(const TayType& x) { 
#ifdef TEUCHOS_DEBUG
	TEUCHOS_TEST_FOR_EXCEPTION(is_tay_real(x) == false, std::runtime_error,
			   "Real component is not a differentiable "
			   "function of complex inputs.");
#endif
	TayType y = x;
	y.copyForWrite();
	y.val() = Teuchos::ScalarTraits<ValueT>::real(x.val());
	return y;
      }

      // Imaginary part is only defined for real derivative components
      static TayType imag(const TayType& x) { 
#ifdef TEUCHOS_DEBUG
	TEUCHOS_TEST_FOR_EXCEPTION(is_tay_real(x) == false, std::runtime_error,
			   "Imaginary component is not a differentiable "
			   "function of complex inputs.");
#endif
	return TayType(Teuchos::ScalarTraits<ValueT>::imag(x.val()));
      }

      static ValueT nan() {
	return Teuchos::ScalarTraits<ValueT>::nan(); 
      }
      static bool isnaninf(const TayType& x) { 
	for (int i=0; i<=x.degree(); i++)
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
	return Sacado::StringName<TayType>::eval(); 
      }
      static TayType squareroot(const TayType& x) {
#ifdef TEUCHOS_DEBUG
	TEUCHOS_SCALAR_TRAITS_NAN_INF_ERR(
	  x, "Error, the input value to squareroot(...) a = " << x << 
	  " can not be NaN!" );
#endif
	return std::sqrt(x); 
      }
      static TayType pow(const TayType& x, const TayType& y) { 
	return std::pow(x,y); 
      }

      // Helper function to determine whether a complex value is real
      static bool is_complex_real(const ValueT& x) {
	return 
	  Teuchos::ScalarTraits<ValueT>::magnitude(x-Teuchos::ScalarTraits<ValueT>::real(x)) == 0;
      }

      // Helper function to determine whether a Fad type is real
      static bool is_tay_real(const TayType& x) {
	if (x.size() == 0)
	  return true;
	if (Teuchos::ScalarTraits<ValueT>::isComplex) {
	  for (int i=0; i<=x.degree(); i++)
	    if (!is_complex_real(x.fastAccessCoeff(i)))
	      return false;
	}
	return true;
      }

    }; // class ScalarTraitsImp

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

#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_FAD_SCALARTRAITSIMP_HPP
