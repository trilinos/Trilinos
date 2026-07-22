// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_HAS_EQUAL_TO_HPP
#define SACADO_MPL_HAS_EQUAL_TO_HPP

#include <type_traits>

namespace Sacado {

  namespace mpl {

    template <typename T1, typename T2 = T1, typename = std::void_t<> >
    struct has_equal_to : std::false_type {};

    template <typename T1, typename T2>
    struct has_equal_to<T1, T2, std::void_t<decltype(std::declval<T1>() ==
                                                     std::declval<T2>())> >
    : std::true_type {};

  }

}

#endif // SACADO_MPL_HAS_EQUAL_TO_HPP
