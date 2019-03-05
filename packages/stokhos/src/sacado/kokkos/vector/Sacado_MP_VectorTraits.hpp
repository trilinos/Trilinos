// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_MP_VECTOR_TRAITS_HPP
#define SACADO_MP_VECTOR_TRAITS_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"

//
// Currently only the type-style traits classes will work on the device
// since the function-style traits don't have the needed specializations
// for the device.  This is in particular true for the scalar traits.
//

// Forward declarations
namespace Sacado {
  namespace MP {
    template <typename S> class Vector;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Vector types
  SACADO_AD_PROMOTE_SPEC( MP, Vector )

  //! Specialization of %ScalarType to Vector types
  template <typename S>
  struct ScalarType< MP::Vector<S> > {
    typedef typename ScalarType<typename MP::Vector<S>::value_type>::type type;
  };

  //! Specialization of %ValueType to Vector types
  template <typename S>
  struct ValueType< MP::Vector<S> > {
    typedef typename MP::Vector<S>::value_type type;
  };

  //! Specialization of %IsADType to Vector types
  template <typename S>
  struct IsADType< MP::Vector<S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Vector types
  template <typename S>
  struct IsScalarType< MP::Vector<S> > {
    static const bool value = S::is_static;
  };

  //! Specialization of %Value to Vector types
  template <typename S>
  struct Value< MP::Vector<S> > {
    typedef typename ValueType< MP::Vector<S> >::type value_type;
    KOKKOS_INLINE_FUNCTION
    static const value_type& eval(const MP::Vector<S>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Vector types
  template <typename S>
  struct ScalarValue< MP::Vector<S> > {
    typedef typename ValueType< MP::Vector<S> >::type value_type;
    typedef typename ScalarType< MP::Vector<S> >::type scalar_type;
    KOKKOS_INLINE_FUNCTION
    static const scalar_type& eval(const MP::Vector<S>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Vector types
  template <typename S>
  struct StringName< MP::Vector<S> > {
    static std::string eval() {
      return std::string("Sacado::MP::Vector< ") +
        StringName<S>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to Vector types
  template <typename S>
  struct IsEqual< MP::Vector<S> > {
    KOKKOS_INLINE_FUNCTION
    static bool eval(const MP::Vector<S>& x,
                     const MP::Vector<S>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to Vector types
  template <typename S>
  struct IsStaticallySized< MP::Vector<S> > {
    static const bool value = S::is_static;
  };

  //! Specialization of %StaticSize to Vector types
  template <typename S>
  struct StaticSize< MP::Vector<S> > {
    static const unsigned value = S::static_size;
  };

} // namespace Sacado

// Define Teuchos traits classes
// Note:  Stokhos has required dependency on all Teuchos sub-packages
#include "Stokhos_ConfigDefs.h"
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_MP_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_as.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to Vector types
  template <typename S>
  struct PromotionTraits< Sacado::MP::Vector<S>,
                          Sacado::MP::Vector<S> > {
    typedef typename Sacado::Promote< Sacado::MP::Vector<S>,
                                      Sacado::MP::Vector<S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to Vector types
  template <typename S, typename R>
  struct PromotionTraits< Sacado::MP::Vector<S>, R > {
    typedef typename Sacado::Promote< Sacado::MP::Vector<S>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to Vector types
  template <typename L, typename S>
  struct PromotionTraits< L, Sacado::MP::Vector<S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::MP::Vector<S> >::type
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
#if defined(HAVE_STOKHOS_ENSEMBLE_REDUCT)
  template <typename S>
  struct ScalarTraits< Sacado::MP::Vector<S> > :
    public Sacado::MP::ScalarTraitsImp<S,true> {};
#else
  template <typename S>
  struct ScalarTraits< Sacado::MP::Vector<S> > :
    public Sacado::MP::ScalarTraitsImp<S,false> {};
#endif

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename S>
  struct SerializationTraits<Ordinal, Sacado::MP::Vector<S> > :
    public Sacado::MP::SerializationTraitsImp< Ordinal,
                                                Sacado::MP::Vector<S>,
                                                S::is_static > {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::MP::Vector<S> > :
    public Sacado::MP::SerializerImp< Ordinal,
                                       Sacado::MP::Vector<S>,
                                       ValueTypeSerializer<Ordinal,typename Sacado::MP::Vector<S>::value_type> >
  {
    typedef Sacado::MP::Vector<S> VecType;
    typedef typename VecType::value_type value_type;
    typedef ValueTypeSerializer<Ordinal,value_type> ValueSerializer;
    typedef Sacado::MP::SerializerImp< Ordinal,VecType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };

//! Specializations for Teuchos::as<T>
template<class TypeTo, class StorageFrom>
class ValueTypeConversionTraits< TypeTo, Sacado::MP::Vector<StorageFrom> > {
public:
  typedef Sacado::MP::Vector<StorageFrom> TypeFrom;
  //! Convert t from a TypeFrom object to a TypeTo object.
  static TypeTo convert (const TypeFrom& t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.
    return Teuchos::as<TypeTo>(t.coeff(0));
  }

  //! Convert t from a TypeFrom object to a TypeTo object, with checks for validity.
  static TypeTo safeConvert (const TypeFrom& t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.  No
    // runtime checking (e.g., for overflow) can be done by default;
    // only specializations can define meaningful and portable
    // run-time checks of conversions.
    return Teuchos::as<TypeTo>(t.coeff(0));
  }
};

template<class TypeTo, class ExprFrom>
class ValueTypeConversionTraits< TypeTo, Sacado::MP::Expr<ExprFrom> > {
public:
  typedef Sacado::MP::Expr<ExprFrom> TypeFrom;
  //! Convert t from a TypeFrom object to a TypeTo object.
  static TypeTo convert (const TypeFrom& t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.
    return Teuchos::as<TypeTo>(t.derived().coeff(0));
  }

  //! Convert t from a TypeFrom object to a TypeTo object, with checks for validity.
  static TypeTo safeConvert (const TypeFrom& t) {
    // This default implementation is just an implicit conversion and
    // may generate compiler warnings on dangerous conversions.  No
    // runtime checking (e.g., for overflow) can be done by default;
    // only specializations can define meaningful and portable
    // run-time checks of conversions.
    return Teuchos::as<TypeTo>(t.derived().coeff(0));
  }
};

// Should also do TypeTo, and TypeTo,TypeFrom as MP::Vector, but the real way
// to fix is to make sure it is never called at all (requires fixing
// magnitudeType)

}

#endif // SACADO_MP_VECTORTRAITS_HPP
