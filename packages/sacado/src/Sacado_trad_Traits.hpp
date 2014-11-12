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

#ifndef SACADO_TRAD_TRAITS_HPP
#define SACADO_TRAD_TRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Rad {
    template <typename T> class ADvar;
    template <typename T> class ADvari;
  }
}

namespace Sacado {

  //! Specialization of %Promote to ADvar types
  SACADO_AD_PROMOTE_SPEC( Rad::ADvar )

  //! Specialization of %ScalarType to ADvar types
  template <typename T>
  struct ScalarType< Rad::ADvar<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ScalarType to ADvari types
  template <typename T>
  struct ScalarType< Rad::ADvari<T> > {
    typedef typename ScalarType<T>::type type;
  };

  //! Specialization of %ValueType to ADvar types
  template <typename T>
  struct ValueType< Rad::ADvar<T> > {
    typedef T type;
  };

  //! Specialization of %ValueType to ADvari types
  template <typename T>
  struct ValueType< Rad::ADvari<T> > {
    typedef T type;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsADType< Rad::ADvar<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsADType< Rad::ADvari<T> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ADvar types
  template <typename T>
  struct IsScalarType< Rad::ADvar<T> > {
    static const bool value = false;
  };

  //! Specialization of %IsADType to ADvari types
  template <typename T>
  struct IsScalarType< Rad::ADvari<T> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ADvar types
  template <typename T>
  struct Value< Rad::ADvar<T> > {
    typedef typename ValueType< Rad::ADvar<T> >::type value_type;
    static value_type eval(const Rad::ADvar<T>& x) {
      return x.val(); }
  };

  //! Specialization of %MarkConstant to ADvar types
  template <typename T>
  struct MarkConstant< Rad::ADvar<T> > {
    static void eval(Rad::ADvar<T>& x) { AD_Const(x); }
  };

  //! Specialization of %MarkConstant to ADvari types
  template <typename T>
  struct MarkConstant< Rad::ADvari<T> > {
    static void eval(Rad::ADvari<T>& x) { AD_Const(x); }
  };

  //! Specialization of %ScalarValue to ADvar types
  template <typename T>
  struct ScalarValue< Rad::ADvar<T> > {
    typedef typename ValueType< Rad::ADvar<T> >::type value_type;
    typedef typename ScalarType< Rad::ADvar<T> >::type scalar_type;
    static scalar_type eval(const Rad::ADvar<T>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to ADvar types
  template <typename T>
  struct StringName< Rad::ADvar<T> > {
    static std::string eval() {
      return std::string("Sacado::Rad::ADvar< ") +
        StringName<T>::eval() + " >"; }
  };

} // namespace Sacado

#endif // SACADO_TRAD_TRAITS_HPP
