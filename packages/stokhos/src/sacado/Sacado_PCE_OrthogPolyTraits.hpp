// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_PCE_ORTHOGPOLYTRAITS_HPP
#define SACADO_PCE_ORTHOGPOLYTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace PCE {
    template <typename T> class OrthogPoly;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T>
  class Promote< PCE::OrthogPoly<T>, PCE::OrthogPoly<T> > {
  public:

    typedef PCE::OrthogPoly<T> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R>
  class Promote< PCE::OrthogPoly<L>, R > {
  public:

    typedef typename ValueType< PCE::OrthogPoly<L> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::OrthogPoly<value_type> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R>
  class Promote< L, PCE::OrthogPoly<R> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< PCE::OrthogPoly<R> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef PCE::OrthogPoly<value_type> type;
  };

  //! Specialization of %ScalarType to OrthogPoly types
  template <typename T>
  struct ScalarType< PCE::OrthogPoly<T> > {
    typedef typename ScalarType<typename PCE::OrthogPoly<T>::value_type>::type type;
  };

  //! Specialization of %ValueType to OrthogPoly types
  template <typename T>
  struct ValueType< PCE::OrthogPoly<T> > {
    typedef typename PCE::OrthogPoly<T>::value_type type;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T>
  struct IsADType< PCE::OrthogPoly<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T>
  struct IsScalarType< PCE::OrthogPoly<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to OrthogPoly types
  template <typename T>
  struct Value< PCE::OrthogPoly<T> > {
    typedef typename ValueType< PCE::OrthogPoly<T> >::type value_type;
    static const value_type& eval(const PCE::OrthogPoly<T>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to OrthogPoly types
  template <typename T>
  struct ScalarValue< PCE::OrthogPoly<T> > {
    typedef typename ValueType< PCE::OrthogPoly<T> >::type value_type;
    typedef typename ScalarType< PCE::OrthogPoly<T> >::type scalar_type;
    static const scalar_type& eval(const PCE::OrthogPoly<T>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to OrthogPoly types
  template <typename T>
  struct StringName< PCE::OrthogPoly<T> > {
    static std::string eval() { 
      return std::string("Sacado::PCE::OrthogPoly< ") + 
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
  template <typename T>
  struct PromotionTraits< Sacado::PCE::OrthogPoly<T>, 
			  Sacado::PCE::OrthogPoly<T> > {
    typedef typename Sacado::Promote< Sacado::PCE::OrthogPoly<T>,
				      Sacado::PCE::OrthogPoly<T> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename R>
  struct PromotionTraits< Sacado::PCE::OrthogPoly<T>, R > {
    typedef typename Sacado::Promote< Sacado::PCE::OrthogPoly<T>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T>
  struct PromotionTraits< L, Sacado::PCE::OrthogPoly<T> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::PCE::OrthogPoly<T> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename T>
  struct ScalarTraits< Sacado::PCE::OrthogPoly<T> > {
    typedef Sacado::PCE::OrthogPoly<T> ScalarType;
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

  }; // class ScalarTraits< Sacado::PCE::OrthogPoly<T> >
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_PCE_UNIVARIATEHERMITETRAITS_HPP
