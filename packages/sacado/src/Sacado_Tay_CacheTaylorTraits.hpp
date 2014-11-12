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

#ifndef SACADO_TAYLOR_DTAYLORTRAITS_HPP
#define SACADO_TAYLOR_DTAYLORTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Tay {
    template <typename T> class CacheTaylor;
  }
}

namespace Sacado {

  //! Specialization of %Promote to CacheTaylor types
  SACADO_AD_PROMOTE_SPEC( Tay::CacheTaylor )

  //! Specialization of %ScalarType to DFad types
  template <typename T>
  struct ScalarType< Tay::CacheTaylor<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to DFad types
  template <typename T>
  struct ValueType< Tay::CacheTaylor<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to DFad types
  template <typename T>
  struct IsADType< Tay::CacheTaylor<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to DFad types
  template <typename T>
  struct IsScalarType< Tay::CacheTaylor<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to DFad types
  template <typename T>
  struct Value< Tay::CacheTaylor<T> > {
    typedef typename ValueType< Tay::CacheTaylor<T> >::type value_type;
    static const value_type& eval(const Tay::CacheTaylor<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to CacheTaylor types
  template <typename T>
  struct ScalarValue< Tay::CacheTaylor<T> > {
    typedef typename ValueType< Tay::CacheTaylor<T> >::type value_type;
    typedef typename ScalarType< Tay::CacheTaylor<T> >::type scalar_type;
    static const scalar_type& eval(const Tay::CacheTaylor<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to CacheTaylor types
  template <typename T>
  struct StringName< Tay::CacheTaylor<T> > {
    static std::string eval() {
      return std::string("Sacado::Tay::CacheTaylor< ") +
        StringName<T>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to Taylor types
  template <typename T>
  struct IsEqual< Tay::CacheTaylor<T> > {
    static bool eval(const Tay::CacheTaylor<T>& x,
                     const Tay::CacheTaylor<T>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to Taylor types
  template <typename T>
  struct IsStaticallySized< Tay::CacheTaylor<T> > {
    static const bool value = false;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Tay_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT>
  struct PromotionTraits< Sacado::Tay::CacheTaylor<ValueT>,
                          Sacado::Tay::CacheTaylor<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Tay::CacheTaylor<ValueT>,
                                      Sacado::Tay::CacheTaylor<ValueT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Tay::CacheTaylor<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Tay::CacheTaylor<ValueT>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Tay::CacheTaylor<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Tay::CacheTaylor<ValueT> >::type
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename ValueT>
  struct ScalarTraits< Sacado::Tay::CacheTaylor<ValueT> > :
    public Sacado::Tay::ScalarTraitsImp< Sacado::Tay::CacheTaylor<ValueT> >
  {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename ValueT>
  struct SerializationTraits<Ordinal, Sacado::Tay::CacheTaylor<ValueT> > :
    public Sacado::Tay::SerializationTraitsImp< Ordinal,
                                                Sacado::Tay::CacheTaylor<ValueT> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename ValueT>
  struct ValueTypeSerializer<Ordinal, Sacado::Tay::CacheTaylor<ValueT> > :
    public Sacado::Tay::SerializerImp< Ordinal,
                                       Sacado::Tay::CacheTaylor<ValueT>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::Tay::CacheTaylor<ValueT> TayType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Tay::SerializerImp< Ordinal,TayType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_TAYLOR_DTAYLORTRAITS_HPP
