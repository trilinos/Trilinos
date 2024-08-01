// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_FOR_EACH_HPP
#define SACADO_MPL_FOR_EACH_HPP

#include "Sacado_ConfigDefs.h"

#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_deref.hpp"

namespace Sacado {

  namespace mpl {

    template <class Seq,
              class Iter1 = typename mpl::begin<Seq>::type,
              class Iter2 = typename mpl::end<Seq>::type>
    struct for_each {
      template <typename Op>
      SACADO_INLINE_FUNCTION
      for_each(const Op& op) {
        op(typename mpl::deref<Iter1>::type());
        for_each<Seq, typename mpl::next<Iter1>::type, Iter2> f(op);
      }
    };

    template <class Seq, class Iter1>
    struct for_each<Seq, Iter1, Iter1> {
      template <typename Op>
      SACADO_INLINE_FUNCTION
      for_each(const Op& op) {}
    };

    // Same as for_each above, but without SACADO_INLINE_FUNCTION for functors
    // that aren't meant to run inside kokkos kernels.
    template <class Seq,
              class Iter1 = typename mpl::begin<Seq>::type,
              class Iter2 = typename mpl::end<Seq>::type>
    struct for_each_no_kokkos {
      template <typename Op>
      for_each_no_kokkos(const Op& op) {
        op(typename mpl::deref<Iter1>::type());
        for_each_no_kokkos<Seq, typename mpl::next<Iter1>::type, Iter2> f(op);
      }
    };

    template <class Seq, class Iter1>
    struct for_each_no_kokkos<Seq, Iter1, Iter1> {
      template <typename Op>
      for_each_no_kokkos(const Op& /* op */) {}
    };

  }

}

#endif // SACADO_MPL_FOR_EACH_HPP
