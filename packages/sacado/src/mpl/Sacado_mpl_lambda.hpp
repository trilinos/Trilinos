// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_LAMBDA_HPP
#define SACADO_MPL_LAMBDA_HPP

#include "Sacado_mpl_bind.hpp"
#include "Sacado_mpl_quote.hpp"
#include "Sacado_mpl_type_wrap.hpp"
#include "Sacado_mpl_if.hpp"
#include "Sacado_mpl_is_placeholder.hpp"

namespace Sacado {

  namespace mpl {

    template <class F> struct lambda : mpl::type_wrap<F> {};

    template <template<class T1> class F,
              class T1>
    struct lambda< F<T1> > :
      mpl_if< is_placeholder< F<T1> >,
              type_wrap< bind1< quote1<F>,
                                typename lambda<T1>::type > >,
              type_wrap< F<T1> > > {};

    template <template<class T1, class T2> class F,
              class T1,
              class T2>
    struct lambda< F<T1,T2> > :
      mpl_if< is_placeholder< F<T1,T2> >,
              type_wrap< bind2< quote2<F>,
                                typename lambda<T1>::type,
                                typename lambda<T2>::type > >,
              type_wrap< F<T1,T2> > > {};

    template <template<class T1, class T2, class T3> class F,
              class T1,
              class T2,
              class T3>
    struct lambda< F<T1,T2,T3> > :
      mpl_if< is_placeholder< F<T1,T2,T3> >,
              type_wrap< bind3< quote3<F>,
                                typename lambda<T1>::type,
                                typename lambda<T2>::type,
                                typename lambda<T3>::type > >,
              type_wrap< F<T1,T2,T3> > > {};

    template <template<class T1, class T2, class T3, class T4> class F,
              class T1,
              class T2,
              class T3,
              class T4>
    struct lambda< F<T1,T2,T3,T4> > :
      mpl_if< is_placeholder< F<T1,T2,T3,T4> >,
              type_wrap< bind4< quote4<F>,
                                typename lambda<T1>::type,
                                typename lambda<T2>::type,
                                typename lambda<T3>::type,
                                typename lambda<T4>::type > >,
              type_wrap< F<T1,T2,T3,T4> > > {};

    template <template<class T1, class T2, class T3, class T4, class T5> class F,
              class T1,
              class T2,
              class T3,
              class T4,
              class T5>
    struct lambda< F<T1,T2,T3,T4,T5> > :
      mpl_if< is_placeholder< F<T1,T2,T3,T4,T5> >,
              type_wrap< bind5< quote5<F>,
                                typename lambda<T1>::type,
                                typename lambda<T2>::type,
                                typename lambda<T3>::type,
                                typename lambda<T4>::type,
                                typename lambda<T5>::type > >,
              type_wrap< F<T1,T2,T3,T4,T5> > > {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_LAMBDA_HPP
