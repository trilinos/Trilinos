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
// @HEADER

#ifndef SACADO_ELRCACHEFAD_SLFADTRAITS_HPP
#define SACADO_ELRCACHEFAD_SLFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace ELRCacheFad {
    template <typename T, int Num> class SLFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SLFad types
  template <typename ValueT, int Num>
  struct Promote< ELRCacheFad::SLFad<ValueT,Num>,
                  ELRCacheFad::SLFad<ValueT,Num> > {
    typedef ELRCacheFad::SLFad<ValueT,Num> type;
  };

  //! Specialization of %Promote to SLFad types
  template <typename ValueT, int Num, typename R>
  struct Promote< ELRCacheFad::SLFad<ValueT,Num>, R > {
    typedef typename ValueType< ELRCacheFad::SLFad<ValueT,Num> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRCacheFad::SLFad<value_type,Num> type;
  };

  //! Specialization of %Promote to SLFad types
  template <typename L, typename ValueT, int Num>
  struct Promote< L, ELRCacheFad::SLFad<ValueT, Num> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ELRCacheFad::SLFad<ValueT,Num> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRCacheFad::SLFad<value_type,Num> type;
  };

  //! Specialization of %ScalarType to SLFad types
  template <typename ValueT, int Num>
  struct ScalarType< ELRCacheFad::SLFad<ValueT,Num> > {
    typedef typename ELRCacheFad::SLFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SLFad types
  template <typename ValueT, int Num>
  struct ValueType< ELRCacheFad::SLFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num>
  struct IsADType< ELRCacheFad::SLFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num>
  struct IsScalarType< ELRCacheFad::SLFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SLFad types
  template <typename ValueT, int Num>
  struct Value< ELRCacheFad::SLFad<ValueT,Num> > {
    typedef typename ValueType< ELRCacheFad::SLFad<ValueT,Num> >::type value_type;
    KOKKOS_INLINE_FUNCTION
    static const value_type& eval(const ELRCacheFad::SLFad<ValueT,Num>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SLFad types
  template <typename ValueT, int Num>
  struct ScalarValue< ELRCacheFad::SLFad<ValueT,Num> > {
    typedef typename ValueType< ELRCacheFad::SLFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< ELRCacheFad::SLFad<ValueT,Num> >::type scalar_type;
    KOKKOS_INLINE_FUNCTION
    static const scalar_type& eval(const ELRCacheFad::SLFad<ValueT,Num>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SLFad types
  template <typename ValueT, int Num>
  struct StringName< ELRCacheFad::SLFad<ValueT,Num> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Sacado::ELRCacheFad::SLFad< "
         << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str();
    }
  };

  //! Specialization of %IsEqual to DFad types
  template <typename ValueT, int Num>
  struct IsEqual< ELRCacheFad::SLFad<ValueT,Num> > {
    KOKKOS_INLINE_FUNCTION
    static bool eval(const ELRCacheFad::SLFad<ValueT,Num>& x,
                     const ELRCacheFad::SLFad<ValueT,Num>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SLFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< ELRCacheFad::SLFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %IsStaticallySized to SLFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< const ELRCacheFad::SLFad<ValueT,Num> > {
    static const bool value = false;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, int Num>
  struct PromotionTraits< Sacado::ELRCacheFad::SLFad<ValueT,Num>,
                          Sacado::ELRCacheFad::SLFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::ELRCacheFad::SLFad<ValueT,Num>,
                                      Sacado::ELRCacheFad::SLFad<ValueT,Num> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::ELRCacheFad::SLFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::ELRCacheFad::SLFad<ValueT,Num>,
                                      R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::ELRCacheFad::SLFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L,
                                      Sacado::ELRCacheFad::SLFad<ValueT,Num> >::type
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::ELRCacheFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::ELRCacheFad::SLFad<ValueT,Num> >
  {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::ELRCacheFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal,
                                                Sacado::ELRCacheFad::SLFad<ValueT,Num> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::ELRCacheFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal,
                                       Sacado::ELRCacheFad::SLFad<ValueT,Num>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::ELRCacheFad::SLFad<ValueT,Num> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ELRCACHEFAD_SLFADTRAITS_HPP
