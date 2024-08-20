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

#ifndef SACADO_CACHEFAD_DFADTRAITS_HPP
#define SACADO_CACHEFAD_DFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace CacheFad {
    template <typename T> class DFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to DFad types
  SACADO_FAD_PROMOTE_SPEC( CacheFad, DFad )

  //! Specialization of %ScalarType to DFad types
  template <typename ValueT>
  struct ScalarType< CacheFad::DFad<ValueT> > {
    typedef typename CacheFad::DFad<ValueT>::ScalarT type;
  };

  //! Specialization of %ValueType to DFad types
  template <typename ValueT>
  struct ValueType< CacheFad::DFad<ValueT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to DFad types
  template <typename ValueT>
  struct IsADType< CacheFad::DFad<ValueT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to DFad types
  template <typename ValueT>
  struct IsScalarType< CacheFad::DFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to DFad types
  template <typename ValueT>
  struct Value< CacheFad::DFad<ValueT> > {
    typedef typename ValueType< CacheFad::DFad<ValueT> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const CacheFad::DFad<ValueT>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DFad types
  template <typename ValueT>
  struct ScalarValue< CacheFad::DFad<ValueT> > {
    typedef typename ValueType< CacheFad::DFad<ValueT> >::type value_type;
    typedef typename ScalarType< CacheFad::DFad<ValueT> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const CacheFad::DFad<ValueT>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DFad types
  template <typename ValueT>
  struct StringName< CacheFad::DFad<ValueT> > {
    static std::string eval() {
      return std::string("Sacado::CacheFad::DFad< ") +
        StringName<ValueT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to DFad types
  template <typename ValueT>
  struct IsEqual< CacheFad::DFad<ValueT> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const CacheFad::DFad<ValueT>& x,
                     const CacheFad::DFad<ValueT>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to DFad types
  template <typename ValueT>
  struct IsStaticallySized< CacheFad::DFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %IsStaticallySized to DFad types
  template <typename ValueT>
  struct IsStaticallySized< const CacheFad::DFad<ValueT> > {
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
  struct PromotionTraits< Sacado::CacheFad::DFad<ValueT>,
                          Sacado::CacheFad::DFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::CacheFad::DFad<ValueT>,
                                      Sacado::CacheFad::DFad<ValueT> >::type
    promote;
  };

  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::CacheFad::DFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::CacheFad::DFad<ValueT>, R >::type
    promote;
  };

  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::CacheFad::DFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::CacheFad::DFad<ValueT> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT>
  struct ScalarTraits< Sacado::CacheFad::DFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::CacheFad::DFad<ValueT> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT>
  struct SerializationTraits<Ordinal, Sacado::CacheFad::DFad<ValueT> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::CacheFad::DFad<ValueT> >
  {};

  template <typename Ordinal, typename ValueT>
  struct ValueTypeSerializer<Ordinal, Sacado::CacheFad::DFad<ValueT> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::CacheFad::DFad<ValueT>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::CacheFad::DFad<ValueT> FadType;
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

#endif // SACADO_FAD_DFADTRAITS_HPP
