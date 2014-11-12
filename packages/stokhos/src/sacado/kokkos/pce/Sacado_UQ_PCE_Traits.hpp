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

#ifndef SACADO_UQ_PCE_TRAITS_HPP
#define SACADO_UQ_PCE_TRAITS_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"

// Forward declarations
namespace Sacado {
  namespace UQ {
    template <typename S> class PCE;
  }
}

namespace Sacado {

  //! Specialization of %Promote to PCE types
  template <typename S>
  struct Promote< UQ::PCE<S>,
                  UQ::PCE<S> > {
    typedef UQ::PCE<S> type;
  };
  template <typename S>
  struct Promote< UQ::PCE<S>,
                  typename UQ::PCE<S>::value_type > {
    typedef UQ::PCE<S> type;
  };
  template <typename S>
  struct Promote< typename UQ::PCE<S>::value_type,
                  UQ::PCE<S> > {
    typedef UQ::PCE<S> type;
  };
  template <typename S>
  struct Promote< UQ::PCE<S>,
                  typename dummy< typename UQ::PCE<S>::value_type,
                                  typename UQ::PCE<S>::scalar_type
                                  >::type > {
    typedef UQ::PCE<S> type;
  };
  template <typename S>
  struct Promote< typename dummy< typename UQ::PCE<S>::value_type,
                                  typename UQ::PCE<S>::scalar_type
                                  >::type,
                  UQ::PCE<S> > {
    typedef UQ::PCE<S> type;
  };

  //! Specialization of %ScalarType to PCE types
  template <typename S>
  struct ScalarType< UQ::PCE<S> > {
    typedef typename ScalarType<typename UQ::PCE<S>::value_type>::type type;
  };

  //! Specialization of %ValueType to PCE types
  template <typename S>
  struct ValueType< UQ::PCE<S> > {
    typedef typename UQ::PCE<S>::value_type type;
  };

  //! Specialization of %IsADType to PCE types
  template <typename S>
  struct IsADType< UQ::PCE<S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to PCE types
  template <typename S>
  struct IsScalarType< UQ::PCE<S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to PCE types
  template <typename S>
  struct Value< UQ::PCE<S> > {
    typedef typename ValueType< UQ::PCE<S> >::type value_type;
    static const value_type& eval(const UQ::PCE<S>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to PCE types
  template <typename S>
  struct ScalarValue< UQ::PCE<S> > {
    typedef typename ValueType< UQ::PCE<S> >::type value_type;
    typedef typename ScalarType< UQ::PCE<S> >::type scalar_type;
    static const scalar_type& eval(const UQ::PCE<S>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to PCE types
  template <typename S>
  struct StringName< UQ::PCE<S> > {
    static std::string eval() {
      return std::string("Sacado::UQ::PCE< ") +
        StringName<S>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to PCE types
  template <typename S>
  struct IsEqual< UQ::PCE<S> > {
    static bool eval(const UQ::PCE<S>& x,
                     const UQ::PCE<S>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to PCE types
  template <typename S>
  struct IsStaticallySized< UQ::PCE<S> > {
    static const bool value = S::is_static;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_UQ_PCE_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename S>
  struct PromotionTraits< Sacado::UQ::PCE<S>,
                          Sacado::UQ::PCE<S> > {
    typedef typename Sacado::Promote< Sacado::UQ::PCE<S>,
                                      Sacado::UQ::PCE<S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename S, typename R>
  struct PromotionTraits< Sacado::UQ::PCE<S>, R > {
    typedef typename Sacado::Promote< Sacado::UQ::PCE<S>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename S>
  struct PromotionTraits< L, Sacado::UQ::PCE<S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::UQ::PCE<S> >::type
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename S>
  struct ScalarTraits< Sacado::UQ::PCE<S> > :
    public Sacado::UQ::PCEScalarTraitsImp< Sacado::UQ::PCE<S> > {};

  //! Specializtion of %Teuchos::ValueTypeConversionTraits
  template <typename TypeTo, typename S>
  struct ValueTypeConversionTraits< TypeTo, Sacado::UQ::PCE<S> > :
    public Sacado::UQ::PCEValueTypeConversionTraitsImp< TypeTo,
                                                        Sacado::UQ::PCE<S> > {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename S>
  struct SerializationTraits<Ordinal, Sacado::UQ::PCE<S> > :
    public Sacado::UQ::PCESerializationTraitsImp< Ordinal,
                                                  Sacado::UQ::PCE<S> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::UQ::PCE<S> > :
    public Sacado::UQ::PCESerializerImp< Ordinal,
                                         Sacado::UQ::PCE<S>,
                                         ValueTypeSerializer<Ordinal,typename S::value_type> >
  {
    typedef Sacado::UQ::PCE<S> PCEType;
    typedef ValueTypeSerializer<Ordinal,typename S::value_type> ValueSerializer;
    typedef Sacado::UQ::PCESerializerImp< Ordinal,PCEType,ValueSerializer> Base;
    typedef typename Base::cijk_type cijk_type;
    ValueTypeSerializer(const cijk_type& cijk,
                        const Teuchos::RCP<const ValueSerializer>& vs) :
      Base(cijk,vs) {}
  };

// Should also do TypeTo, and TypeTo,TypeFrom as UQ::PCE, but the real way
// to fix is to make sure it is never called at all (requires fixing
// magnitudeType)

}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_UQ_PCE_TRAITS_HPP
