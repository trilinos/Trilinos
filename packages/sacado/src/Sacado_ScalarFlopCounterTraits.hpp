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
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
// @HEADER

#ifndef SACADO_SCALARFLOPCOUNTERTRAITS_HPP
#define SACADO_SCALARFLOPCOUNTERTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace FlopCounterPack {
    template <typename T> class ScalarFlopCounter;
  }
}

namespace Sacado {

  //! Specialization of %Promote to ScalarFlopCounter types
  SACADO_AD_PROMOTE_SPEC( FlopCounterPack, ScalarFlopCounter )

  //! Specialization of %ScalarType to ScalarFlopCounter types
  template <typename ScalarT>
  struct ScalarType< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    typedef typename ScalarType<ScalarT>::type type;
  };

  //! Specialization of %ValueType to ScalarFlopCounter types
  template <typename ScalarT>
  struct ValueType< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    typedef ScalarT type;
  };

  //! Specialization of %IsADType to ScalarFlopCounter types
  template <typename ScalarT>
  struct IsADType< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ScalarFlopCounter types
  template <typename ScalarT>
  struct IsScalarType< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ScalarFlopCounter types
  template <typename ScalarT>
  struct Value< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    typedef typename ValueType< FlopCounterPack::ScalarFlopCounter<ScalarT> >::type value_type;
    static const value_type& eval(const FlopCounterPack::ScalarFlopCounter<ScalarT>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to ScalarFlopCounter types
  template <typename ScalarT>
  struct ScalarValue< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    typedef typename ValueType< FlopCounterPack::ScalarFlopCounter<ScalarT> >::type value_type;
    typedef typename ScalarType< FlopCounterPack::ScalarFlopCounter<ScalarT> >::type scalar_type;
    static const scalar_type& eval(const FlopCounterPack::ScalarFlopCounter<ScalarT>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to ScalarFlopCounter types
  template <typename ScalarT>
  struct StringName< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    static std::string eval() {
      return std::string("Sacado::FlopCounterPack::ScalarFlopCounter< ") +
        StringName<ScalarT>::eval() + " >"; }
  };

  // Need specialization for IsEqual

  //! Specialization of %IsStaticallySized to ScalarFlopCounter types
  template <typename ScalarT>
  struct IsStaticallySized< FlopCounterPack::ScalarFlopCounter<ScalarT> > {
    static const bool value = true;
  };

} // namespace Sacado

#endif // SACADO_SCALARFLOPCOUNTERTRAITS_HPP
