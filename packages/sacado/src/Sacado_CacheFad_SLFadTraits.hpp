// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_SLFADTRAITS_HPP
#define SACADO_CACHEFAD_SLFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace CacheFad {
    template <typename T, int Num> class SLFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SLFad types
  SACADO_SFAD_PROMOTE_SPEC( CacheFad, SLFad )

  //! Specialization of %ScalarType to SLFad types
  template <typename ValueT, int Num>
  struct ScalarType< CacheFad::SLFad<ValueT,Num> > {
    typedef typename CacheFad::SLFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SLFad types
  template <typename ValueT, int Num>
  struct ValueType< CacheFad::SLFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num>
  struct IsADType< CacheFad::SLFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num>
  struct IsScalarType< CacheFad::SLFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SLFad types
  template <typename ValueT, int Num>
  struct Value< CacheFad::SLFad<ValueT,Num> > {
    typedef typename ValueType< CacheFad::SLFad<ValueT,Num> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const CacheFad::SLFad<ValueT,Num>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SLFad types
  template <typename ValueT, int Num>
  struct ScalarValue< CacheFad::SLFad<ValueT,Num> > {
    typedef typename ValueType< CacheFad::SLFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< CacheFad::SLFad<ValueT,Num> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const CacheFad::SLFad<ValueT,Num>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SLFad types
  template <typename ValueT, int Num>
  struct StringName< CacheFad::SLFad<ValueT,Num> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Sacado::CacheFad::SLFad< "
         << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str();
    }
  };

  //! Specialization of %IsEqual to DFad types
  template <typename ValueT, int Num>
  struct IsEqual< CacheFad::SLFad<ValueT,Num> > {
     SACADO_INLINE_FUNCTION
    static bool eval(const CacheFad::SLFad<ValueT,Num>& x,
                     const CacheFad::SLFad<ValueT,Num>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SLFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< CacheFad::SLFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %IsStaticallySized to SLFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< const CacheFad::SLFad<ValueT,Num> > {
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
  template <typename ValueT, int Num>
  struct PromotionTraits< Sacado::CacheFad::SLFad<ValueT,Num>,
                          Sacado::CacheFad::SLFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::CacheFad::SLFad<ValueT,Num>,
                                      Sacado::CacheFad::SLFad<ValueT,Num> >::type
    promote;
  };

  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::CacheFad::SLFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::CacheFad::SLFad<ValueT,Num>, R >::type
    promote;
  };

  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::CacheFad::SLFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::CacheFad::SLFad<ValueT,Num> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::CacheFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::CacheFad::SLFad<ValueT,Num> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::CacheFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::CacheFad::SLFad<ValueT,Num> >
  {};

  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::CacheFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::CacheFad::SLFad<ValueT,Num>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::CacheFad::SLFad<ValueT,Num> FadType;
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

#endif // SACADO_CACHEFAD_SLFADTRAITS_HPP
