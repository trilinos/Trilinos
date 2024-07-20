// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_GENERALFADTRAITS_HPP
#define SACADO_FAD_EXP_GENERALFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
  namespace Exp {
    template <typename S> class GeneralFad;
  }
  }
}

namespace Sacado {

  //! Specialization of %Promote to GeneralFad types
  SACADO_FAD_PROMOTE_SPEC( Fad::Exp, GeneralFad )

  //! Specialization of %ScalarType to GeneralFad types
  template <typename Storage>
  struct ScalarType< Fad::Exp::GeneralFad<Storage> > {
    typedef typename Fad::Exp::GeneralFad<Storage>::scalar_type type;
  };

  //! Specialization of %Storageype to GeneralFad types
  template <typename Storage>
  struct ValueType< Fad::Exp::GeneralFad<Storage> > {
    typedef typename Fad::Exp::GeneralFad<Storage>::value_type type;
  };

  //! Specialization of %IsADType to GeneralFad types
  template <typename Storage>
  struct IsADType< Fad::Exp::GeneralFad<Storage> > {
    static const bool value = true;
  };

  //! Specialization of %IsScalarType to GeneralFad types
  template <typename Storage>
  struct IsScalarType< Fad::Exp::GeneralFad<Storage> > {
    static const bool value = false;
  };

  //! Specialization of %IsSimdType to GeneralFad types
  template <typename Storage>
  struct IsSimdType< Fad::Exp::GeneralFad<Storage> > {
    static const bool value =
      IsSimdType< typename Fad::Exp::GeneralFad<Storage>::value_type >::value;
  };

  //! Specialization of %Value to GeneralFad types
  template <typename Storage>
  struct Value< Fad::Exp::GeneralFad<Storage> > {
    typedef typename ValueType< Fad::Exp::GeneralFad<Storage> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const Fad::Exp::GeneralFad<Storage>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to GeneralFad types
  template <typename Storage>
  struct ScalarValue< Fad::Exp::GeneralFad<Storage> > {
    typedef typename ValueType< Fad::Exp::GeneralFad<Storage> >::type value_type;
    typedef typename ScalarType< Fad::Exp::GeneralFad<Storage> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const Fad::Exp::GeneralFad<Storage>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to GeneralFad types
  template <typename Storage>
  struct StringName< Fad::Exp::GeneralFad<Storage> > {
    static std::string eval() {
      return std::string("Sacado::Fad::Exp::GeneralFad< ") +
        StringName<typename Storage::value_type>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to GeneralFad types
  template <typename Storage>
  struct IsEqual< Fad::Exp::GeneralFad<Storage> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const Fad::Exp::GeneralFad<Storage>& x,
                     const Fad::Exp::GeneralFad<Storage>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to GeneralFad types
  template <typename Storage>
  struct IsStaticallySized< Fad::Exp::GeneralFad<Storage> > {
    static const bool value = Storage::is_statically_sized;
  };

  //! Specialization of %IsStaticallySized to GeneralFad types
  template <typename Storage>
  struct IsStaticallySized< const Fad::Exp::GeneralFad<Storage> > {
    static const bool value = Storage::is_statically_sized;
  };

  //! Specialization of %StaticSize to GeneralFad types
  template <typename Storage>
  struct StaticSize< Fad::Exp::GeneralFad<Storage> > {
    static const unsigned value = Storage::static_size;
  };

  //! Specialization of %StaticSize to GeneralFad types
  template <typename Storage>
  struct StaticSize< const Fad::Exp::GeneralFad<Storage> > {
    static const unsigned value = Storage::static_size;
  };

} // namespace Sacado

//
// Define Teuchos traits classes
//

// Promotion traits
#ifdef HAVE_SACADO_TEUCHOSNUMERICS
#include "Teuchos_PromotionTraits.hpp"
namespace Teuchos {
  template <typename Storage>
  struct PromotionTraits< Sacado::Fad::Exp::GeneralFad<Storage>,
                          Sacado::Fad::Exp::GeneralFad<Storage> > {
    typedef typename Sacado::Promote< Sacado::Fad::Exp::GeneralFad<Storage>,
                                      Sacado::Fad::Exp::GeneralFad<Storage> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to GeneralFad types
  template <typename Storage, typename R>
  struct PromotionTraits< Sacado::Fad::Exp::GeneralFad<Storage>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::Exp::GeneralFad<Storage>,
                                      R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to GeneralFad types
  template <typename L, typename Storage>
  struct PromotionTraits< L, Sacado::Fad::Exp::GeneralFad<Storage> > {
  public:
    typedef typename Sacado::Promote< L,
                                      Sacado::Fad::Exp::GeneralFad<Storage> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename Storage>
  struct ScalarTraits< Sacado::Fad::Exp::GeneralFad<Storage> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::Exp::GeneralFad<Storage> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename Storage>
  struct SerializationTraits<Ordinal, Sacado::Fad::Exp::GeneralFad<Storage> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::Fad::Exp::GeneralFad<Storage> >
  {};

  template <typename Ordinal, typename Storage>
  struct ValueTypeSerializer<Ordinal, Sacado::Fad::Exp::GeneralFad<Storage> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::Fad::Exp::GeneralFad<Storage>,
                                       ValueTypeSerializer<Ordinal,typename Storage::value_type> >
  {
    typedef Sacado::Fad::Exp::GeneralFad<Storage> FadType;
    typedef ValueTypeSerializer<Ordinal,typename Storage::value_type> ValueSerializer;
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

#endif // SACADO_FAD_EXP_GENERALFADTRAITS_HPP
