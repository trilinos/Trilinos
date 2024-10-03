// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TRADVEC_TRAITS_HPP
#define SACADO_TRADVEC_TRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace RadVec {
    template <typename T> class ADvar;
    template <typename T> class ADvari;
  }
}

namespace Sacado {

  //! Specialization of %Promote to ADvar types
  SACADO_RAD_PROMOTE_SPEC( RadVec )

  //! Specialization of %ScalarType to ADvar types
  template <typename T>
  struct ScalarType< RadVec::ADvar<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ScalarType to ADvari types
  template <typename T>
  struct ScalarType< RadVec::ADvari<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to ADvar types
  template <typename T>
  struct ValueType< RadVec::ADvar<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to ADvari types
  template <typename T>
  struct ValueType< RadVec::ADvari<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsADType< RadVec::ADvar<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsADType< RadVec::ADvari<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsScalarType< RadVec::ADvar<T> > {
    static const bool value = false;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsScalarType< RadVec::ADvari<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ADvar types
  template <typename T>
  struct Value< RadVec::ADvar<T> > {
    typedef typename ValueType< RadVec::ADvar<T> >::type value_type;
    static value_type eval(const RadVec::ADvar<T>& x) {
      return x.val(); }
  };

  //! Specialization of %MarkConstant to ADvar types
  template <typename T>
  struct MarkConstant< RadVec::ADvar<T> > {
    static void eval(RadVec::ADvar<T>& x) { AD_Const(x); }
  };

  //! Specialization of %MarkConstant to ADvari types
  template <typename T>
  struct MarkConstant< RadVec::ADvari<T> > {
    static void eval(RadVec::ADvari<T>& x) { AD_Const(x); }
  };

  //! Specialization of %ScalarValue to ADvar types
  template <typename T>
  struct ScalarValue< RadVec::ADvar<T> > {
    typedef typename ValueType< RadVec::ADvar<T> >::type value_type;
    typedef typename ScalarType< RadVec::ADvar<T> >::type scalar_type;
    static scalar_type eval(const RadVec::ADvar<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to ADvar types
  template <typename T>
  struct StringName< RadVec::ADvar<T> > {
    static std::string eval() {
      return std::string("Sacado::RadVec::ADvar< ") +
        StringName<T>::eval() + " >"; }
  };

} // namespace Sacado

#endif // SACADO_TRADVEC_TRAITS_HPP
