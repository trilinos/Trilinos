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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.  
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
// @HEADER

#ifndef SACADO_ELRFAD_DFADTRAITS_HPP
#define SACADO_ELRFAD_DFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ELRFad {
    template <typename T1, typename T2> class DFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to DFad types
  template <typename ValueT, typename ScalarT>
  struct Promote< ELRFad::DFad<ValueT,ScalarT>, ELRFad::DFad<ValueT,ScalarT> > {
    typedef ELRFad::DFad<ValueT,ScalarT> type;
  };

  //! Specialization of %Promote to DFad types
  template <typename ValueT, typename ScalarT, typename R>
  struct Promote< ELRFad::DFad<ValueT,ScalarT>, R > {
    typedef typename ValueType< ELRFad::DFad<ValueT,ScalarT> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRFad::DFad<value_type,ScalarT> type;
  };

  //! Specialization of %Promote to DFad types
  template <typename L, typename ValueT, typename ScalarT>
  struct Promote< L, ELRFad::DFad<ValueT, ScalarT> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ELRFad::DFad<ValueT,ScalarT> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRFad::DFad<value_type,ScalarT> type;
  };

  //! Specialization of %ScalarType to DFad types
  template <typename ValueT, typename ScalarT>
  struct ScalarType< ELRFad::DFad<ValueT,ScalarT> > {
    typedef ScalarT type;
  };

  //! Specialization of %ValueType to DFad types
  template <typename ValueT, typename ScalarT>
  struct ValueType< ELRFad::DFad<ValueT,ScalarT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to DFad types
  template <typename ValueT, typename ScalarT>
  struct IsADType< ELRFad::DFad<ValueT,ScalarT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to DFad types
  template <typename ValueT, typename ScalarT>
  struct IsScalarType< ELRFad::DFad<ValueT,ScalarT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to DFad types
  template <typename ValueT, typename ScalarT>
  struct Value< ELRFad::DFad<ValueT,ScalarT> > {
    typedef typename ValueType< ELRFad::DFad<ValueT,ScalarT> >::type value_type;
    static const value_type& eval(const ELRFad::DFad<ValueT,ScalarT>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DELRFad types
  template <typename ValueT, typename ScalarT>
  struct ScalarValue< ELRFad::DFad<ValueT,ScalarT> > {
    typedef typename ValueType< ELRFad::DFad<ValueT,ScalarT> >::type value_type;
    typedef typename ScalarType< ELRFad::DFad<ValueT,ScalarT> >::type scalar_type;
    static const scalar_type& eval(const ELRFad::DFad<ValueT,ScalarT>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DFad types
  template <typename ValueT, typename ScalarT>
  struct StringName< ELRFad::DFad<ValueT,ScalarT> > {
    static std::string eval() { 
      return std::string("Sacado::ELRFad::DFad< ") + 
	StringName<ValueT>::eval() + ", " + 
	StringName<ScalarT>::eval() + " >"; }
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, typename ScalarT>
  struct PromotionTraits< Sacado::ELRFad::DFad<ValueT,ScalarT>, 
			  Sacado::ELRFad::DFad<ValueT,ScalarT> > {
    typedef typename Sacado::Promote< Sacado::ELRFad::DFad<ValueT,ScalarT>,
				      Sacado::ELRFad::DFad<ValueT,ScalarT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, typename ScalarT, typename R>
  struct PromotionTraits< Sacado::ELRFad::DFad<ValueT,ScalarT>, R > {
    typedef typename Sacado::Promote< Sacado::ELRFad::DFad<ValueT,ScalarT>,
				      R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename ValueT, typename ScalarT>
  struct PromotionTraits< L, Sacado::ELRFad::DFad<ValueT, ScalarT> > {
  public:
    typedef typename Sacado::Promote< L, 
				      Sacado::ELRFad::DFad<ValueT,ScalarT> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename ValueT, typename ScalarT>
  struct ScalarTraits< Sacado::ELRFad::DFad<ValueT,ScalarT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::ELRFad::DFad<ValueT,ScalarT> >
  {};
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ELRFAD_DFADTRAITS_HPP
