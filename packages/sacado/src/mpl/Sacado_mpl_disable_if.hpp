// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_DISABLE_IF_HPP
#define SACADO_MPL_DISABLE_IF_HPP

namespace Sacado {

  namespace mpl {

    template <bool, typename T = void>
    struct disable_if_c {};

    template <typename T>
    struct disable_if_c<false, T> {
      typedef T type;
    };

    template <class Cond, typename T = void>
    struct disable_if
      : disable_if_c<Cond::value, T> {};

    template <bool, typename T = void>
    struct lazy_disable_if_c {};

    template <typename T>
    struct lazy_disable_if_c<false, T> {
      typedef typename T::type type;
    };

    template <class Cond, typename T = void>
    struct lazy_disable_if
      : lazy_disable_if_c<Cond::value, T> {};

  }

}

#endif // SACADO_MPL_DISABLE_IF_HPP
