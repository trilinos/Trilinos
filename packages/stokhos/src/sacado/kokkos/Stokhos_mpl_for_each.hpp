// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MPL_FOR_EACH_HPP
#define STOKHOS_MPL_FOR_EACH_HPP

#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_deref.hpp"

#include "Kokkos_Macros.hpp"

namespace Stokhos {

  namespace mpl {

    template <class Seq,
              class Iter1 = typename Sacado::mpl::begin<Seq>::type,
              class Iter2 = typename Sacado::mpl::end<Seq>::type>
    struct for_each {
      template <typename Op>
      KOKKOS_INLINE_FUNCTION
      for_each(const Op& op) {
        op(typename Sacado::mpl::deref<Iter1>::type());
        for_each<Seq, typename Sacado::mpl::next<Iter1>::type, Iter2> f(op);
      }
    };

    template <class Seq, class Iter1>
    struct for_each<Seq, Iter1, Iter1> {
      template <typename Op>
      KOKKOS_INLINE_FUNCTION
      for_each(const Op& op) {}
    };

  }

}

#endif // STOKHOS_MPL_FOR_EACH_HPP
