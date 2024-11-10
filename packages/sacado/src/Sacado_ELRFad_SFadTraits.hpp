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

#ifndef SACADO_ELRFAD_SFADTRAITS_HPP
#define SACADO_ELRFAD_SFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace ELRFad {
    template <typename T, int Num> class SFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SFad types
  SACADO_SFAD_PROMOTE_SPEC( ELRFad, SFad )

  //! Specialization of %ScalarType to SFad types
  template <typename ValueT, int Num>
  struct ScalarType< ELRFad::SFad<ValueT,Num> > {
    typedef typename ELRFad::SFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SFad types
  template <typename ValueT, int Num>
  struct ValueType< ELRFad::SFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsADType< ELRFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsScalarType< ELRFad::SFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SFad types
  template <typename ValueT, int Num>
  struct Value< ELRFad::SFad<ValueT,Num> > {
    typedef typename ValueType< ELRFad::SFad<ValueT,Num> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const ELRFad::SFad<ValueT,Num>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SFad types
  template <typename ValueT, int Num>
  struct ScalarValue< ELRFad::SFad<ValueT,Num> > {
    typedef typename ValueType< ELRFad::SFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< ELRFad::SFad<ValueT,Num> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const ELRFad::SFad<ValueT,Num>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SFad types
  template <typename ValueT, int Num>
  struct StringName< ELRFad::SFad<ValueT,Num> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Sacado::ELRFad::SFad< "
         << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str();
    }
  };

  //! Specialization of %IsEqual to SFad types
  template <typename ValueT, int Num>
  struct IsEqual< ELRFad::SFad<ValueT,Num> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const ELRFad::SFad<ValueT,Num>& x,
                     const ELRFad::SFad<ValueT,Num>& y) {
       return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< ELRFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< const ELRFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< ELRFad::SFad<ValueT,Num> > {
    static const unsigned value = Num;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< const ELRFad::SFad<ValueT,Num> > {
    static const unsigned value = Num;
  };

} // namespace Sacado

//
// Define Teuchos traits classes
//

// Promotion traits
#ifdef HAVE_SACADO_TEUCHOSNUMERICS
#include "Teuchos_PromotionTraits.hpp"
namespace Teuchos {
  template <typename ValueT, int Num>
  struct PromotionTraits< Sacado::ELRFad::SFad<ValueT,Num>,
                          Sacado::ELRFad::SFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::ELRFad::SFad<ValueT,Num>,
                                      Sacado::ELRFad::SFad<ValueT,Num> >::type
    promote;
  };

  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::ELRFad::SFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::ELRFad::SFad<ValueT,Num>, R >::type
    promote;
  };

  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::ELRFad::SFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ELRFad::SFad<ValueT,Num> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::ELRFad::SFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::ELRFad::SFad<ValueT,Num> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::ELRFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::ELRFad::SFad<ValueT,Num> >
  {};

  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::ELRFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::ELRFad::SFad<ValueT,Num>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::ELRFad::SFad<ValueT,Num> FadType;
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

#endif // SACADO_ELRFAD_SFADTRAITS_HPP
