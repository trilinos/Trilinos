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

#ifndef SACADO_ELRCACHEFAD_SFADTRAITS_HPP
#define SACADO_ELRCACHEFAD_SFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace ELRCacheFad {
    template <typename T, int Num> class SFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SFad types
  SACADO_SFAD_PROMOTE_SPEC( ELRCacheFad, SFad )

  //! Specialization of %ScalarType to SFad types
  template <typename ValueT, int Num>
  struct ScalarType< ELRCacheFad::SFad<ValueT,Num> > {
    typedef typename ELRCacheFad::SFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SFad types
  template <typename ValueT, int Num>
  struct ValueType< ELRCacheFad::SFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsADType< ELRCacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsScalarType< ELRCacheFad::SFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SFad types
  template <typename ValueT, int Num>
  struct Value< ELRCacheFad::SFad<ValueT,Num> > {
    typedef typename ValueType< ELRCacheFad::SFad<ValueT,Num> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const ELRCacheFad::SFad<ValueT,Num>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SFad types
  template <typename ValueT, int Num>
  struct ScalarValue< ELRCacheFad::SFad<ValueT,Num> > {
    typedef typename ValueType< ELRCacheFad::SFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< ELRCacheFad::SFad<ValueT,Num> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const ELRCacheFad::SFad<ValueT,Num>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SFad types
  template <typename ValueT, int Num>
  struct StringName< ELRCacheFad::SFad<ValueT,Num> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Sacado::ELRCacheFad::SFad< "
         << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str();
    }
  };

  //! Specialization of %IsEqual to SFad types
  template <typename ValueT, int Num>
  struct IsEqual< ELRCacheFad::SFad<ValueT,Num> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const ELRCacheFad::SFad<ValueT,Num>& x,
                     const ELRCacheFad::SFad<ValueT,Num>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< ELRCacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< const ELRCacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< ELRCacheFad::SFad<ValueT,Num> > {
    static const unsigned value = Num;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< const ELRCacheFad::SFad<ValueT,Num> > {
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
  struct PromotionTraits< Sacado::ELRCacheFad::SFad<ValueT,Num>,
                          Sacado::ELRCacheFad::SFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::ELRCacheFad::SFad<ValueT,Num>,
                                      Sacado::ELRCacheFad::SFad<ValueT,Num> >::type
    promote;
  };

  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::ELRCacheFad::SFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::ELRCacheFad::SFad<ValueT,Num>, R >::type
    promote;
  };

  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::ELRCacheFad::SFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ELRCacheFad::SFad<ValueT,Num> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Fad_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::ELRCacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::ELRCacheFad::SFad<ValueT,Num> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Fad_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::ELRCacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::ELRCacheFad::SFad<ValueT,Num> >
  {};

  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::ELRCacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::ELRCacheFad::SFad<ValueT,Num>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::ELRCacheFad::SFad<ValueT,Num> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif

// KokkosComm
#if defined(HAVE_SACADO_KOKKOSCORE) && defined(HAVE_SACADO_TEUCHOSKOKKOSCOMM) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
#include "KokkosExp_View_Fad.hpp"
#endif

#endif // SACADO_ELRCACHEFAD_SFADTRAITS_HPP
