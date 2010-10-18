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

#ifndef SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
#define SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETPCE {
    template <typename T, typename S> class OrthogPoly;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T, typename S>
  class Promote< ETPCE::OrthogPoly<T,S>, ETPCE::OrthogPoly<T,S> > {
  public:

    typedef ETPCE::OrthogPoly<T,S> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R, typename S>
  class Promote< ETPCE::OrthogPoly<L,S>, R > {
  public:

    typedef typename ValueType< ETPCE::OrthogPoly<L,S> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETPCE::OrthogPoly<value_type,S> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R, typename S>
  class Promote< L, ETPCE::OrthogPoly<R,S> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ETPCE::OrthogPoly<R,S> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETPCE::OrthogPoly<value_type,S> type;
  };

  //! Specialization of %ScalarType to OrthogPoly types
  template <typename T, typename S>
  struct ScalarType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ScalarType<typename ETPCE::OrthogPoly<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to OrthogPoly types
  template <typename T, typename S>
  struct ValueType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ETPCE::OrthogPoly<T,S>::value_type type;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsADType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsScalarType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to OrthogPoly types
  template <typename T, typename S>
  struct Value< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    static const value_type& eval(const ETPCE::OrthogPoly<T,S>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to OrthogPoly types
  template <typename T, typename S>
  struct ScalarValue< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    typedef typename ScalarType< ETPCE::OrthogPoly<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETPCE::OrthogPoly<T,S>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to OrthogPoly types
  template <typename T, typename S>
  struct StringName< ETPCE::OrthogPoly<T,S> > {
    static std::string eval() { 
      return std::string("Sacado::ETPCE::OrthogPoly< ") + 
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
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, 
			  Sacado::ETPCE::OrthogPoly<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>,
				      Sacado::ETPCE::OrthogPoly<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETPCE::OrthogPoly<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETPCE::OrthogPoly<T,S> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETPCE::OrthogPoly<T,S> > {
    typedef Sacado::ETPCE::OrthogPoly<T,S> ScalarType;
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
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(is_pce_real(a), std::runtime_error,
			 "Complex conjugate is not defined for "
			 "complex PCE inputs.");
#endif
      //return std::fabs(a); 
      magnitudeType b(a.size());
      if (Teuchos::ScalarTraits<ValueT>::real(a.val()) >= 0)
	for (int i=0; i<a.size(); i++)
	  b.fastAccessCoeff(i) = 
	    Teuchos::ScalarTraits<ValueT>::magnitude(a.fastAccessCoeff(i));
      else
	for (int i=0; i<a.size(); i++)
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
    static ScalarType conjugate(const ScalarType& x) {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(is_pce_real(x), std::runtime_error,
			 "Complex conjugate is not defined for "
			 "complex PCE inputs.");
#endif
	return x;
    }   
    
    // Real part is only defined for real derivative components
    static ScalarType real(const ScalarType& x) { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(is_pce_real(x) == false, std::runtime_error,
			 "Real component is not defined for "
			 "complex PCE inputs.");
#endif
      return x;
    }
    
    // Imaginary part is only defined for real derivative components
    static ScalarType imag(const ScalarType& x) { 
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(is_pce_real(x) == false, std::runtime_error,
			 "Imaginary component is not defined for "
			 "complex PCE inputs.");
#endif
      return x;
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

    // Helper function to determine whether a complex value is real
    static bool is_complex_real(const ValueT& x) {
      return 
	Teuchos::ScalarTraits<ValueT>::magnitude(x-Teuchos::ScalarTraits<ValueT>::real(x)) == 0;
    }

    // Helper function to determine whether a Fad type is real
    static bool is_pce_real(const ScalarType& x) {
      if (x.size() == 0)
	return true;
      if (Teuchos::ScalarTraits<ValueT>::isComplex) {
	for (int i=0; i<x.size(); i++)
	  if (!is_complex_real(x.fastAccessCoeff(i)))
	    return false;
      }
      return true;
    }

  }; // class ScalarTraits< Sacado::ETPCE::OrthogPoly<T,S> >
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
