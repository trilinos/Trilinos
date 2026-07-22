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

#ifndef SACADO_FAD_SIMPLEFADTRAITS_HPP
#define SACADO_FAD_SIMPLEFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T> class SimpleFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SimpleFad types
  SACADO_FAD_PROMOTE_SPEC( Fad, SimpleFad )

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

  //! Specialization of %IsScalarType to SimpleFad types
  template <typename ValueT>
  struct IsScalarType< Fad::SimpleFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %IsSimdType to SimpleFad types
  template <typename ValueT>
  struct IsSimdType< Fad::SimpleFad<ValueT> > {
    static const bool value = IsSimdType<ValueT>::value;
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

  //! Specialization of %IsEqual to SimpleFad types
  template <typename ValueT>
  struct IsEqual< Fad::SimpleFad<ValueT> > {
    static bool eval(const Fad::SimpleFad<ValueT>& x,
                     const Fad::SimpleFad<ValueT>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SimpleFad types
  template <typename ValueT>
  struct IsStaticallySized< Fad::SimpleFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %IsFad to SimpleFad types
  template <typename ValueT>
  struct IsFad< Fad::SimpleFad<ValueT> > {
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
  struct PromotionTraits< Sacado::Fad::SimpleFad<ValueT>,
                          Sacado::Fad::SimpleFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Fad::SimpleFad<ValueT>,
                                      Sacado::Fad::SimpleFad<ValueT> >::type
    promote;
  };

  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Fad::SimpleFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::SimpleFad<ValueT>, R >::type
    promote;
  };

  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Fad::SimpleFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::SimpleFad<ValueT> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT>
  struct ScalarTraits< Sacado::Fad::SimpleFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::SimpleFad<ValueT> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT>
  struct SerializationTraits<Ordinal, Sacado::Fad::SimpleFad<ValueT> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::Fad::SimpleFad<ValueT> >
  {};

  template <typename Ordinal, typename ValueT>
  struct ValueTypeSerializer<Ordinal, Sacado::Fad::SimpleFad<ValueT> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::Fad::SimpleFad<ValueT>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::Fad::SimpleFad<ValueT> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif

#endif // SACADO_FAD_SIMPLEFADTRAITS_HPP
