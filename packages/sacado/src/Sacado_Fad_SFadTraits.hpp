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

#ifndef SACADO_FAD_SFADTRAITS_HPP
#define SACADO_FAD_SFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T, int Num> class SFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SFad types
  template <typename ValueT, int Num>
  struct Promote< Fad::SFad<ValueT,Num>,
                  Fad::SFad<ValueT,Num> > {
    typedef Fad::SFad<ValueT,Num> type;
  };

  //! Specialization of %Promote to SFad types
  template <typename ValueT, int Num, typename R>
  struct Promote< Fad::SFad<ValueT,Num>, R > {
    typedef typename ValueType< Fad::SFad<ValueT,Num> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::SFad<value_type,Num> type;
  };

  //! Specialization of %Promote to SFad types
  template <typename L, typename ValueT, int Num>
  struct Promote< L, Fad::SFad<ValueT, Num> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< Fad::SFad<ValueT,Num> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::SFad<value_type,Num> type;
  };

  //! Specialization of %ScalarType to SFad types
  template <typename ValueT, int Num>
  struct ScalarType< Fad::SFad<ValueT,Num> > {
    typedef typename Fad::SFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SFad types
  template <typename ValueT, int Num>
  struct ValueType< Fad::SFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsADType< Fad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsScalarType< Fad::SFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SFad types
  template <typename ValueT, int Num>
  struct Value< Fad::SFad<ValueT,Num> > {
    typedef typename ValueType< Fad::SFad<ValueT,Num> >::type value_type;
    KOKKOS_INLINE_FUNCTION
    static const value_type& eval(const Fad::SFad<ValueT,Num>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SFad types
  template <typename ValueT, int Num>
  struct ScalarValue< Fad::SFad<ValueT,Num> > {
    typedef typename ValueType< Fad::SFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< Fad::SFad<ValueT,Num> >::type scalar_type;
    KOKKOS_INLINE_FUNCTION
    static const scalar_type& eval(const Fad::SFad<ValueT,Num>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SFad types
  template <typename ValueT, int Num>
  struct StringName< Fad::SFad<ValueT,Num> > {
    static std::string eval() {
       std::stringstream ss;
      ss << "Sacado::Fad::SFad< "
         << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str();
    }
  };

  //! Specialization of %IsEqual to SFad types
  template <typename ValueT, int Num>
  struct IsEqual< Fad::SFad<ValueT,Num> > {
    KOKKOS_INLINE_FUNCTION
    static bool eval(const Fad::SFad<ValueT,Num>& x,
                     const Fad::SFad<ValueT,Num>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< Fad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< const Fad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< Fad::SFad<ValueT,Num> > {
    static const unsigned value = Num;
  };

  //! Specialization of %StaticSize to SFad types
  template <typename ValueT, int Num>
  struct StaticSize< const Fad::SFad<ValueT,Num> > {
    static const unsigned value = Num;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to SFad types
  template <typename ValueT, int Num>
  struct PromotionTraits< Sacado::Fad::SFad<ValueT,Num>,
                          Sacado::Fad::SFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::Fad::SFad<ValueT,Num>,
                                      Sacado::Fad::SFad<ValueT,Num> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to SFad types
  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::Fad::SFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::SFad<ValueT,Num>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to SFad types
  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::Fad::SFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::SFad<ValueT,Num> >::type
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::Fad::SFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::SFad<ValueT,Num> >
  {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::Fad::SFad<ValueT,Num> > :
    public Sacado::Fad::StaticSerializationTraitsImp< Ordinal,
                                                      Sacado::Fad::SFad<ValueT,Num> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::Fad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::Fad::SFad<ValueT,Num>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::Fad::SFad<ValueT,Num> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_FAD_SFADTRAITS_HPP
