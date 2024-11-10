// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_TRAD2_TRAITS_HPP
#define SACADO_TRAD2_TRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Rad2 {
    template <typename T> class ADvar;
    template <typename T> class ADvari;
  }
}

namespace Sacado {

  //! Specialization of %Promote to ADvar types
  SACADO_RAD_PROMOTE_SPEC( Rad2 )

  //! Specialization of %ScalarType to ADvar types
  template <typename T>
  struct ScalarType< Rad2::ADvar<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ScalarType to ADvari types
  template <typename T>
  struct ScalarType< Rad2::ADvari<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to ADvar types
  template <typename T>
  struct ValueType< Rad2::ADvar<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to ADvari types
  template <typename T>
  struct ValueType< Rad2::ADvari<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsADType< Rad2::ADvar<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsADType< Rad2::ADvari<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsScalarType< Rad2::ADvar<T> > {
    static const bool value = false;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsScalarType< Rad2::ADvari<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ADvar types
  template <typename T>
  struct Value< Rad2::ADvar<T> > {
    typedef typename ValueType< Rad2::ADvar<T> >::type value_type;
    static value_type eval(const Rad2::ADvar<T>& x) {
      return x.val(); }
  };

  //! Specialization of %MarkConstant to ADvar types
  template <typename T>
  struct MarkConstant< Rad2::ADvar<T> > {
    static void eval(Rad2::ADvar<T>& x) { AD_Const(x); }
  };

  //! Specialization of %MarkConstant to ADvari types
  template <typename T>
  struct MarkConstant< Rad2::ADvari<T> > {
    static void eval(Rad2::ADvari<T>& x) { AD_Const(x); }
  };

  //! Specialization of %ScalarValue to ADvar types
  template <typename T>
  struct ScalarValue< Rad2::ADvar<T> > {
    typedef typename ValueType< Rad2::ADvar<T> >::type value_type;
    typedef typename ScalarType< Rad2::ADvar<T> >::type scalar_type;
    static scalar_type eval(const Rad2::ADvar<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to ADvar types
  template <typename T>
  struct StringName< Rad2::ADvar<T> > {
    static std::string eval() {
      return std::string("Sacado::Rad2::ADvar< ") +
        StringName<T>::eval() + " >"; }
  };

} // namespace Sacado

#endif // SACADO_TRAD_TRAITS_HPP
