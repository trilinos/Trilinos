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

#ifndef SACADO_FAD_DVFADTRAITS_HPP
#define SACADO_FAD_DVFADTRAITS_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_GeneralFadTraits.hpp"

#else

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T> class DVFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to DVFad types
  SACADO_FAD_PROMOTE_SPEC( Fad, DVFad )

  //! Specialization of %ScalarType to DVFad types
  template <typename ValueT>
  struct ScalarType< Fad::DVFad<ValueT> > {
    typedef typename Fad::DVFad<ValueT>::ScalarT type;
  };

  //! Specialization of %ValueType to DVFad types
  template <typename ValueT>
  struct ValueType< Fad::DVFad<ValueT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to DVFad types
  template <typename ValueT>
  struct IsADType< Fad::DVFad<ValueT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to DVFad types
  template <typename ValueT>
  struct IsScalarType< Fad::DVFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to DVFad types
  template <typename ValueT>
  struct Value< Fad::DVFad<ValueT> > {
    typedef typename ValueType< Fad::DVFad<ValueT> >::type value_type;
    static const value_type& eval(const Fad::DVFad<ValueT>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DVFad types
  template <typename ValueT>
  struct ScalarValue< Fad::DVFad<ValueT> > {
    typedef typename ValueType< Fad::DVFad<ValueT> >::type value_type;
    typedef typename ScalarType< Fad::DVFad<ValueT> >::type scalar_type;
    static const scalar_type& eval(const Fad::DVFad<ValueT>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DVFad types
  template <typename ValueT>
  struct StringName< Fad::DVFad<ValueT> > {
    static std::string eval() {
      return std::string("Sacado::Fad::DVFad< ") +
        StringName<ValueT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to DVFad types
  template <typename ValueT>
  struct IsEqual< Fad::DVFad<ValueT> > {
    static bool eval(const Fad::DVFad<ValueT>& x, const Fad::DVFad<ValueT>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to DVFad types
  template <typename ValueT>
  struct IsStaticallySized< Fad::DVFad<ValueT> > {
    static const bool value = false;
  };

  template <typename T>
  struct IsFad< Fad::DVFad<T> > {
    static const bool value = true;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DVFad types
  template <typename ValueT>
  struct PromotionTraits< Sacado::Fad::DVFad<ValueT>,
                          Sacado::Fad::DVFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Fad::DVFad<ValueT>,
                                      Sacado::Fad::DVFad<ValueT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DVFad types
  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Fad::DVFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::DVFad<ValueT>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DVFad types
  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Fad::DVFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::DVFad<ValueT> >::type
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename ValueT>
  struct ScalarTraits< Sacado::Fad::DVFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::DVFad<ValueT> >
  {};
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_DVFADTRAITS_HPP
