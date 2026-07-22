// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_TYPE_WRAP_HPP
#define SACADO_MPL_TYPE_WRAP_HPP

#include "Sacado_mpl_has_type.hpp"

namespace Sacado {

  namespace mpl {

    // Wrap a type T so that it becomes a metafunction if it isn't already

    // We should check if T is a class type and derive add_type from T if
    // it is.  This would allow us to use type_wrap in mpl_if since we assume
    // mpl_if<cond,T1,T2> is derived from T1 or T2 in mpl_find.
    template <class T> struct add_type { typedef T type; };

    // Don't use mpl_if because it uses type_wrap in its implementation
    // (actually not, see above).
    template <bool cond, class T> struct type_wrap_impl {};
    template <class T> struct type_wrap_impl<true,T> : T {};
    template <class T> struct type_wrap_impl<false,T> : add_type<T> {};
    template <class T> struct type_wrap :
      type_wrap_impl< mpl::has_type<T>::value, T > {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_TYPE_WRAP_HPP
