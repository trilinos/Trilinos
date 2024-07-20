// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TAY_TAYLORTRAITS_HPP
#define SACADO_TAY_TAYLORTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Tay {
    template <typename T> class Taylor;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  SACADO_AD_PROMOTE_SPEC( Tay, Taylor )

  //! Specialization of %ScalarType to Taylor types
  template <typename T>
  struct ScalarType< Tay::Taylor<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to Taylor types
  template <typename T>
  struct ValueType< Tay::Taylor<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to Taylor types
  template <typename T>
  struct IsADType< Tay::Taylor<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Taylor types
  template <typename T>
  struct IsScalarType< Tay::Taylor<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Taylor types
  template <typename T>
  struct Value< Tay::Taylor<T> > {
    typedef typename ValueType< Tay::Taylor<T> >::type value_type;
    static const value_type& eval(const Tay::Taylor<T>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Taylor types
  template <typename T>
  struct ScalarValue< Tay::Taylor<T> > {
    typedef typename ValueType< Tay::Taylor<T> >::type value_type;
    typedef typename ScalarType< Tay::Taylor<T> >::type scalar_type;
    static const scalar_type& eval(const Tay::Taylor<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Taylor types
  template <typename T>
  struct StringName< Tay::Taylor<T> > {
    static std::string eval() {
      return std::string("Sacado::Tay::Taylor< ") +
        StringName<T>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to Taylor types
  template <typename T>
  struct IsEqual< Tay::Taylor<T> > {
    static bool eval(const Tay::Taylor<T>& x, const Tay::Taylor<T>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to Taylor types
  template <typename T>
  struct IsStaticallySized< Tay::Taylor<T> > {
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
  template <typename ValueT>
  struct PromotionTraits< Sacado::Tay::Taylor<ValueT>,
                          Sacado::Tay::Taylor<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Tay::Taylor<ValueT>,
                                      Sacado::Tay::Taylor<ValueT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Tay::Taylor<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Tay::Taylor<ValueT>, R >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Tay::Taylor<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Tay::Taylor<ValueT> >::type
    promote;
  };
}
#endif

// Scalar traits
#ifdef HAVE_SACADO_TEUCHOSCORE
#include "Sacado_Tay_ScalarTraitsImp.hpp"
namespace Teuchos {
  template <typename ValueT>
  struct ScalarTraits< Sacado::Tay::Taylor<ValueT> > :
    public Sacado::Tay::ScalarTraitsImp< Sacado::Tay::Taylor<ValueT> >
  {};
}
#endif

// Serialization traits
#ifdef HAVE_SACADO_TEUCHOSCOMM
#include "Sacado_Tay_SerializationTraitsImp.hpp"
namespace Teuchos {
  template <typename Ordinal, typename ValueT>
  struct SerializationTraits<Ordinal, Sacado::Tay::Taylor<ValueT> > :
    public Sacado::Tay::SerializationTraitsImp< Ordinal,
                                                Sacado::Tay::Taylor<ValueT> >
  {};

  template <typename Ordinal, typename ValueT>
  struct ValueTypeSerializer<Ordinal, Sacado::Tay::Taylor<ValueT> > :
    public Sacado::Tay::SerializerImp< Ordinal,
                                       Sacado::Tay::Taylor<ValueT>,
                                       ValueTypeSerializer<Ordinal,ValueT> >
  {
    typedef Sacado::Tay::Taylor<ValueT> TayType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Tay::SerializerImp< Ordinal,TayType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
                        Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif

#endif // SACADO_TAYLOR_SIMPLETAYLORTRAITS_HPP
