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

#ifndef SACADO_FAD_SIMPLEFADTRAITS_HPP
#define SACADO_FAD_SIMPLETRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T> class SimpleFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SimpleFad types
  template <typename ValueT>
  struct Promote< Fad::SimpleFad<ValueT>, Fad::SimpleFad<ValueT> > {
    typedef Fad::SimpleFad<ValueT> type;
  };

  //! Specialization of %Promote to SimpleFad types
  template <typename ValueT, typename R>
  struct Promote< Fad::SimpleFad<ValueT>, R > {
    typedef typename ValueType< Fad::SimpleFad<ValueT> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::SimpleFad<value_type> type;
  };

  //! Specialization of %Promote to SimpleFad types
  template <typename L, typename ValueT>
  struct Promote< L, Fad::SimpleFad<ValueT> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< Fad::SimpleFad<ValueT> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::SimpleFad<value_type> type;
  };

  //! Specialization of %ScalarType to SimpleFad types
  template <typename ValueT>
  struct ScalarType< Fad::SimpleFad<ValueT> > {
    typedef typename Fad::SimpleFad<ValueT>::ScalarT type;
  };

  //! Specialization of %ValueType to SimpleFad types
  template <typename ValueT>
  struct ValueType< Fad::SimpleFad<ValueT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SimpleFad types
  template <typename ValueT>
  struct IsADType< Fad::SimpleFad<ValueT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SimpleFad types
  template <typename ValueT>
  struct IsScalarType< Fad::SimpleFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SimpleFad types
  template <typename ValueT>
  struct Value< Fad::SimpleFad<ValueT> > {
    typedef typename ValueType< Fad::SimpleFad<ValueT> >::type value_type;
    static const value_type& eval(const Fad::SimpleFad<ValueT>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SimpleFad types
  template <typename ValueT>
  struct ScalarValue< Fad::SimpleFad<ValueT> > {
    typedef typename ValueType< Fad::SimpleFad<ValueT> >::type value_type;
    typedef typename ScalarType< Fad::SimpleFad<ValueT> >::type scalar_type;
    static const scalar_type& eval(const Fad::SimpleFad<ValueT>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SimpleFad types
  template <typename ValueT>
  struct StringName< Fad::SimpleFad<ValueT> > {
    static std::string eval() { 
      return std::string("Sacado::Fad::SimpleFad< ") + 
	StringName<ValueT>::eval() + " >"; }
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to SimpleFad types
  template <typename ValueT>
  struct PromotionTraits< Sacado::Fad::SimpleFad<ValueT>, 
			  Sacado::Fad::SimpleFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Fad::SimpleFad<ValueT>,
				      Sacado::Fad::SimpleFad<ValueT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to SimpleFad types
  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Fad::SimpleFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::SimpleFad<ValueT>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to SimpleFad types
  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Fad::SimpleFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::SimpleFad<ValueT> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename ValueT>
  struct ScalarTraits< Sacado::Fad::SimpleFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::SimpleFad<ValueT> >
  {};
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_FAD_SIMPLEFADTRAITS_HPP
