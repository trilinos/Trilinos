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

#ifndef SACADO_FAD_DFADTRAITS_HPP
#define SACADO_FAD_DFADTRAITS_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#include "Sacado_Fad_Exp_GeneralFadTraits.hpp"

#else

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T> class DFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to DFad types
  SACADO_FAD_PROMOTE_SPEC( Fad, DFad )

  //! Specialization of %ScalarType to DFad types
  template <typename ValueT>
  struct ScalarType< Fad::DFad<ValueT> > {
    typedef typename Fad::DFad<ValueT>::ScalarT type;
  };

  //! Specialization of %ValueType to DFad types
  template <typename ValueT>
  struct ValueType< Fad::DFad<ValueT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to DFad types
  template <typename ValueT>
  struct IsADType< Fad::DFad<ValueT> > {
    static const bool value = true;
  };

  //! Specialization of %IsScalarType to DFad types
  template <typename ValueT>
  struct IsScalarType< Fad::DFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %IsSimdType to DFad types
  template <typename ValueT>
  struct IsSimdType< Fad::DFad<ValueT> > {
    static const bool value = IsSimdType<ValueT>::value;
  };

  //! Specialization of %Value to DFad types
  template <typename ValueT>
  struct Value< Fad::DFad<ValueT> > {
    typedef typename ValueType< Fad::DFad<ValueT> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const Fad::DFad<ValueT>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DFad types
  template <typename ValueT>
  struct ScalarValue< Fad::DFad<ValueT> > {
    typedef typename ValueType< Fad::DFad<ValueT> >::type value_type;
    typedef typename ScalarType< Fad::DFad<ValueT> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const Fad::DFad<ValueT>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DFad types
  template <typename ValueT>
  struct StringName< Fad::DFad<ValueT> > {
    static std::string eval() {
      return std::string("Sacado::Fad::DFad< ") +
        StringName<ValueT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to DFad types
  template <typename ValueT>
  struct IsEqual< Fad::DFad<ValueT> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const Fad::DFad<ValueT>& x, const Fad::DFad<ValueT>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to DFad types
  template <typename ValueT>
  struct IsStaticallySized< Fad::DFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %IsStaticallySized to DFad types
  template <typename ValueT>
  struct IsStaticallySized< const Fad::DFad<ValueT> > {
    static const bool value = false;
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
  struct PromotionTraits< Sacado::Fad::DFad<ValueT>,
                          Sacado::Fad::DFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Fad::DFad<ValueT>,
                                      Sacado::Fad::DFad<ValueT> >::type
    promote;
  };

  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Fad::DFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::DFad<ValueT>, R >::type
    promote;
  };

  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Fad::DFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::DFad<ValueT> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT>
  struct ScalarTraits< Sacado::Fad::DFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::DFad<ValueT> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT>
  struct SerializationTraits<Ordinal, Sacado::Fad::DFad<ValueT> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::Fad::DFad<ValueT> >
  {};

  template <typename Ordinal, typename ValueT>
  struct ValueTypeSerializer<Ordinal, Sacado::Fad::DFad<ValueT> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::Fad::DFad<ValueT>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::Fad::DFad<ValueT> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif

// KokkosComm
#if defined(HAVE_SACADO_KOKKOS) && defined(HAVE_SACADO_TEUCHOSKOKKOSCOMM) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
#include "KokkosExp_View_Fad.hpp"
#endif

#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // SACADO_FAD_DFADTRAITS_HPP
