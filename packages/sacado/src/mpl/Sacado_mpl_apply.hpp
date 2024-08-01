// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_APPLY_HPP
#define SACADO_MPL_APPLY_HPP

#include "Sacado_mpl_apply_wrap.hpp"
#include "Sacado_mpl_lambda.hpp"
#include "Sacado_mpl_none.hpp"

namespace Sacado {

  namespace mpl {

    template <class F> 
    struct apply0 : apply_wrap0<typename lambda<F>::type> {};

    template <class F, class A1> 
    struct apply1 : apply_wrap1<typename lambda<F>::type,A1> {};

    template <class F, class A1, class A2> 
    struct apply2 : apply_wrap2<typename lambda<F>::type,A1,A2> {};

    template <class F, class A1, class A2, class A3> 
    struct apply3 : apply_wrap3<typename lambda<F>::type,A1,A2,A3> {};

    template <class F, class A1, class A2, class A3, class A4> 
    struct apply4 : apply_wrap4<typename lambda<F>::type,A1,A2,A3,A4> {};

    template <class F, class A1, class A2, class A3, class A4, class A5> 
    struct apply5 : apply_wrap5<typename lambda<F>::type,A1,A2,A3,A4,A5> {};

    template <class F, 
              class A1=mpl::none, 
              class A2=mpl::none, 
              class A3=mpl::none, 
              class A4=mpl::none, 
              class A5=mpl::none>
    struct apply : apply_wrap<typename lambda<F>::type,A1,A2,A3,A4,A5> {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_APPLY_HPP
