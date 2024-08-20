// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
#define SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETPCE {
    template <typename T, typename S> class OrthogPoly;
  }
}

namespace Sacado {

  //! Specialization of %Promote to OrthogPoly types
  SACADO_AD_PROMOTE_SPEC2( ETPCE, OrthogPoly )

  //! Specialization of %ScalarType to OrthogPoly types
  template <typename T, typename S>
  struct ScalarType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ScalarType<typename ETPCE::OrthogPoly<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to OrthogPoly types
  template <typename T, typename S>
  struct ValueType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ETPCE::OrthogPoly<T,S>::value_type type;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsADType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsScalarType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to OrthogPoly types
  template <typename T, typename S>
  struct Value< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    static const value_type& eval(const ETPCE::OrthogPoly<T,S>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to OrthogPoly types
  template <typename T, typename S>
  struct ScalarValue< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    typedef typename ScalarType< ETPCE::OrthogPoly<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETPCE::OrthogPoly<T,S>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to OrthogPoly types
  template <typename T, typename S>
  struct StringName< ETPCE::OrthogPoly<T,S> > {
    static std::string eval() {
      return std::string("Sacado::ETPCE::OrthogPoly< ") +
        StringName<T>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to OrthogPoly types
  template <typename T, typename S>
  struct IsEqual< ETPCE::OrthogPoly<T,S> > {
    static bool eval(const ETPCE::OrthogPoly<T,S>& x,
                     const ETPCE::OrthogPoly<T,S>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to OrthogPoly types
  template <typename T, typename S>
  struct IsStaticallySized< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

} // namespace Sacado

// Define Teuchos traits classes
// Note Stokhos has required dependency on all Teuchos sub-packages
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_PCE_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>,
                          Sacado::ETPCE::OrthogPoly<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>,
                                      Sacado::ETPCE::OrthogPoly<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETPCE::OrthogPoly<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETPCE::OrthogPoly<T,S> >::type
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::ScalarTraitsImp< Sacado::ETPCE::OrthogPoly<T,S> > {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename T, typename S>
  struct SerializationTraits<Ordinal, Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::SerializationTraitsImp< Ordinal,
                                                Sacado::ETPCE::OrthogPoly<T,S> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename T, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::SerializerImp< Ordinal,
                                       Sacado::ETPCE::OrthogPoly<T,S>,
                                       ValueTypeSerializer<Ordinal,T> >
  {
    typedef Sacado::ETPCE::OrthogPoly<T,S> PCEType;
    typedef ValueTypeSerializer<Ordinal,T> ValueSerializer;
    typedef Sacado::PCE::SerializerImp< Ordinal,PCEType,ValueSerializer> Base;
    typedef typename Base::expansion_type expansion_type;
    ValueTypeSerializer(const Teuchos::RCP<expansion_type>& expansion,
                        const Teuchos::RCP<const ValueSerializer>& vs) :
      Base(expansion,vs) {}
  };
}

#endif // SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
