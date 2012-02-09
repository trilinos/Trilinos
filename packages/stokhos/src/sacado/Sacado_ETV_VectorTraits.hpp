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

#ifndef SACADO_ETV_VECTORTRAITS_HPP
#define SACADO_ETV_VECTORTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETV {
    template <typename T, typename S> class Vector;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T, typename S>
  class Promote< ETV::Vector<T,S>, ETV::Vector<T,S> > {
  public:

    typedef ETV::Vector<T,S> type;
  };

  //! Specialization of %Promote to Vector types
  template <typename L, typename R, typename S>
  class Promote< ETV::Vector<L,S>, R > {
  public:

    typedef typename ValueType< ETV::Vector<L,S> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETV::Vector<value_type,S> type;
  };

  //! Specialization of %Promote to Vector types
  template <typename L, typename R, typename S>
  class Promote< L, ETV::Vector<R,S> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ETV::Vector<R,S> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETV::Vector<value_type,S> type;
  };

  //! Specialization of %ScalarType to Vector types
  template <typename T, typename S>
  struct ScalarType< ETV::Vector<T,S> > {
    typedef typename ScalarType<typename ETV::Vector<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to Vector types
  template <typename T, typename S>
  struct ValueType< ETV::Vector<T,S> > {
    typedef typename ETV::Vector<T,S>::value_type type;
  };

  //! Specialization of %IsADType to Vector types
  template <typename T, typename S>
  struct IsADType< ETV::Vector<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Vector types
  template <typename T, typename S>
  struct IsScalarType< ETV::Vector<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Vector types
  template <typename T, typename S>
  struct Value< ETV::Vector<T,S> > {
    typedef typename ValueType< ETV::Vector<T,S> >::type value_type;
    static const value_type& eval(const ETV::Vector<T,S>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Vector types
  template <typename T, typename S>
  struct ScalarValue< ETV::Vector<T,S> > {
    typedef typename ValueType< ETV::Vector<T,S> >::type value_type;
    typedef typename ScalarType< ETV::Vector<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETV::Vector<T,S>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Vector types
  template <typename T, typename S>
  struct StringName< ETV::Vector<T,S> > {
    static std::string eval() { 
      return std::string("Sacado::ETV::Vector< ") + 
	StringName<T>::eval() + " >"; }
  };

  //! Specialization of IsEqual to Vector types
  template <typename T, typename S>
  struct IsEqual< ETV::Vector<T,S> > {
    static bool eval(const ETV::Vector<T,S>& x, 
		     const ETV::Vector<T,S>& y) {
      return x.isEqualTo(y);
    }
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Assert.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S>
  struct PromotionTraits< Sacado::ETV::Vector<T,S>, 
			  Sacado::ETV::Vector<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETV::Vector<T,S>,
				      Sacado::ETV::Vector<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETV::Vector<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETV::Vector<T,S>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETV::Vector<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETV::Vector<T,S> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETV::Vector<T,S> > {
    typedef Sacado::ETV::Vector<T,S> ScalarType;
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

  }; // class ScalarTraits< Sacado::ETV::Vector<T,S> >
}

namespace Sacado {
  namespace ETV {

   //! Serialization implementation for all Vector types
    template <typename Ordinal, typename T, typename S, typename Serializer>
    struct SerializationImp {

    private:

      //! VectorType
      typedef Sacado::ETV::Vector<T,S> VecType;

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

  }

}

namespace Teuchos {

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename T, typename S>
  struct SerializationTraits<Ordinal, Sacado::ETV::Vector<T,S> > {

  private:

    //! VectorType
    typedef Sacado::ETV::Vector<T,S> VecType;

    //! Value type of Vector type
    typedef typename Sacado::ValueType<VecType>::type ValueT;
    
    //! Default serializer for values
    typedef Teuchos::DefaultSerializer<Ordinal,ValueT> DS;

    //! Default serializer type for values
    typedef typename DS::DefaultSerializerType ValueSerializer;

    //! Implementation
    typedef Sacado::ETV::SerializationImp<Ordinal,T,S,ValueSerializer> Imp;

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

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename T, typename S>
  class ValueTypeSerializer<Ordinal, Sacado::ETV::Vector<T,S> > {

  private:

    //! VectorType
    typedef Sacado::ETV::Vector<T,S> VecType;

    //! Value type of Vector type
    typedef typename Sacado::ValueType<VecType>::type ValueT;

    //! Serializer for values
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;

    //! Implementation
    typedef Sacado::ETV::SerializationImp<Ordinal,T,S,ValueSerializer> Imp;
    
    //! Serializer for value types
    Teuchos::RCP<const ValueSerializer> vs;

    //! Specified number of vector components;
    Ordinal sz;
    
  public:

    //! Typename of value serializer
    typedef ValueSerializer value_serializer_type;
    
    /// \brief Whether we support direct serialization.
    static const bool supportsDirectSerialization = 
      Imp::supportsDirectSerialization;
    
    //! Constructor
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs_,
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
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETV_VECTORTRAITS_HPP
