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

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Sacado_mpl_apply.hpp"

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
    typedef typename Sacado::ValueType<T>::type ValueT;
    
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
      ScalarType y(sz);
      for (int i=0; i<sz; i++)
	y.fastAccessCoeff(i) = 
	  Teuchos::ScalarTraits<ValueT>::conjugate(x.fastAccessCoeff(i));
      return y;
    }   
    
    // Real part is only defined for real derivative components
    static ScalarType real(const ScalarType& x) { 
      int sz = x.size();
      ScalarType y(sz);
      for (int i=0; i<sz; i++)
	y.fastAccessCoeff(i) = 
	  Teuchos::ScalarTraits<ValueT>::real(x.fastAccessCoeff(i));
      return y;
    }
    
    // Imaginary part is only defined for real derivative components
    static ScalarType imag(const ScalarType& x) { 
     int sz = x.size();
      ScalarType y(sz);
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
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETV_ORTHOGPOLYTRAITS_HPP
