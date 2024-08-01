// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SCADO_MPL_VECTOR_PUSH_BACK_SPEC_HPP
#define SCADO_MPL_VECTOR_PUSH_BACK_SPEC_HPP

namespace Sacado {

  namespace mpl {

    template <class Vector, class T>
    struct vector_push_back {};

    template <typename...Args, class T>
    struct vector_push_back<mpl::vector<Args...>,T> : mpl::vector<Args...,T> {};

  }

}

#endif // SCADO_MPL_VECTOR_PUSH_BACK_SPEC_HPP
