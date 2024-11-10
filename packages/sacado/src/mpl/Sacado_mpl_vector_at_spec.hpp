// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SCADO_MPL_VECTOR_AT_SPEC_HPP
#define SCADO_MPL_VECTOR_AT_SPEC_HPP

namespace Sacado {

  namespace mpl {

    template <class Vector, int Pos> struct vector_at {};

    template <typename T, typename...Args>
    struct vector_at<mpl::vector<T,Args...>, 0> {
      typedef T type;
    };

    template <typename T, typename...Args, int Pos>
    struct vector_at<mpl::vector<T,Args...>, Pos> {
      typedef typename vector_at<mpl::vector<Args...>, Pos-1>::type type;
    };

  }

}

#endif // SCADO_MPL_VECTOR_AT_SPEC_HPP
