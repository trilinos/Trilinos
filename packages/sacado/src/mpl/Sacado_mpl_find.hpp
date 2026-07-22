// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_FIND_HPP
#define SACADO_MPL_FIND_HPP

#include <type_traits>

#include "Sacado_mpl_none.hpp"
#include "Sacado_mpl_begin.hpp"
#include "Sacado_mpl_end.hpp"
#include "Sacado_mpl_deref.hpp"
#include "Sacado_mpl_next.hpp"
#include "Sacado_mpl_if.hpp"

namespace Sacado {

  namespace mpl {

    template <class Seq, class T>
    class TypeSequenceDoesNotContainType {};

    template <class Seq,
              class T,
              class Iter1 = typename mpl::begin<Seq>::type,
              class Iter2 = typename mpl::end<Seq>::type>
    struct find {
      typedef typename
        mpl::mpl_if< std::is_same<typename mpl::deref<Iter1>::type, T>,
                     Iter1,
                     find<Seq, T, typename mpl::next<Iter1>::type,
                          Iter2> >::type type;
      static const int value = type::value;
    };

    template <class Seq, class T, class Iter1>
    struct find<Seq, T, Iter1, Iter1> {
      static const int value = TypeSequenceDoesNotContainType<Seq,T>::value;
    };

  }

}

#endif // SACADO_MPL_FIND_HPP
