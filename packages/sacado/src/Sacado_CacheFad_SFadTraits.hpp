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
//  A short implementation ( not all operators and
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_SFADTRAITS_HPP
#define SACADO_CACHEFAD_SFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace CacheFad {
    template <typename T, int Num> class SFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SFad types
  SACADO_SFAD_PROMOTE_SPEC( CacheFad, SFad )

  //! Specialization of %ScalarType to SFad types
  template <typename ValueT, int Num>
  struct ScalarType< CacheFad::SFad<ValueT,Num> > {
    typedef typename CacheFad::SFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SFad types
  template <typename ValueT, int Num>
  struct ValueType< CacheFad::SFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsADType< CacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsScalarType< CacheFad::SFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SFad types
  template <typename ValueT, int Num>
  struct Value< CacheFad::SFad<ValueT,Num> > {
    typedef typename ValueType< CacheFad::SFad<ValueT,Num> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const CacheFad::SFad<ValueT,Num>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SFad types
  template <typename ValueT, int Num>
  struct ScalarValue< CacheFad::SFad<ValueT,Num> > {
    typedef typename ValueType< CacheFad::SFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< CacheFad::SFad<ValueT,Num> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const CacheFad::SFad<ValueT,Num>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SFad types
  template <typename ValueT, int Num>
  struct StringName< CacheFad::SFad<ValueT,Num> > {
    static std::string eval() {
       std::stringstream ss;
      ss << "Sacado::CacheFad::SFad< "
         << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str();
    }
  };

  //! Specialization of %IsEqual to SFad types
  template <typename ValueT, int Num>
  struct IsEqual< CacheFad::SFad<ValueT,Num> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const CacheFad::SFad<ValueT,Num>& x,
                     const CacheFad::SFad<ValueT,Num>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< CacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< const CacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< CacheFad::SFad<ValueT,Num> > {
    static const unsigned value = Num;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< const CacheFad::SFad<ValueT,Num> > {
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
  struct PromotionTraits< Sacado::CacheFad::SFad<ValueT,Num>,
                          Sacado::CacheFad::SFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::CacheFad::SFad<ValueT,Num>,
                                      Sacado::CacheFad::SFad<ValueT,Num> >::type
    promote;
  };

  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::CacheFad::SFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::CacheFad::SFad<ValueT,Num>, R >::type
    promote;
  };

  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::CacheFad::SFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::CacheFad::SFad<ValueT,Num> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::CacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::CacheFad::SFad<ValueT,Num> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::CacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::CacheFad::SFad<ValueT,Num> >
  {};

  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::CacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::CacheFad::SFad<ValueT,Num>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::CacheFad::SFad<ValueT,Num> FadType;
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

#endif // SACADO_CACHEFAD_SFADTRAITS_HPP
