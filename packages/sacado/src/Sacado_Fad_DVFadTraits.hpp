// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
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

  //! Specialization of %IsScalarType to DVFad types
  template <typename ValueT>
  struct IsScalarType< Fad::DVFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %IsSimdType to DVFad types
  template <typename ValueT>
  struct IsSimdType< Fad::DVFad<ValueT> > {
    static const bool value = IsSimdType<ValueT>::value;
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

//
// Define Teuchos traits classes
//

// Promotion traits
#ifdef HAVE_SACADO_TEUCHOSNUMERICS
#include "Teuchos_PromotionTraits.hpp"
namespace Teuchos {
  template <typename ValueT>
  struct PromotionTraits< Sacado::Fad::DVFad<ValueT>,
                          Sacado::Fad::DVFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Fad::DVFad<ValueT>,
                                      Sacado::Fad::DVFad<ValueT> >::type
    promote;
  };

  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Fad::DVFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::DVFad<ValueT>, R >::type
    promote;
  };

  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Fad::DVFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::DVFad<ValueT> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT>
  struct ScalarTraits< Sacado::Fad::DVFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::DVFad<ValueT> >
  {};
}
#endif

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_DVFADTRAITS_HPP
