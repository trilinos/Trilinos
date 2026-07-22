// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_QUOTE_HPP
#define SACADO_MPL_QUOTE_HPP

#include "Sacado_mpl_type_wrap.hpp"

namespace Sacado {

  namespace mpl {
   
    // Transform a class/template to a metafunction class 
    // (nested template apply)
    template <class F>
    struct quote0 {
      struct apply : mpl::type_wrap< F > {};
    };

    template <template<class T1> class F>
    struct quote1 {
      template <class T1> 
      struct apply : mpl::type_wrap< F<T1> > {};
    };

    template < template<class T1, 
                        class T2> class F>
    struct quote2 {
      template <class T1, 
                class T2> 
      struct apply : mpl::type_wrap< F<T1,T2> >{};
    };

    template < template<class T1, 
                        class T2, 
                        class T3> class F>
    struct quote3 {
      template <class T1, 
                class T2, 
                class T3> 
      struct apply : mpl::type_wrap< F<T1,T2,T3> >{};
    };

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4> class F>
    struct quote4 {
      template <class T1, 
                class T2, 
                class T3,
                class T4> 
      struct apply : mpl::type_wrap< F<T1,T2,T3,T4> >{};
    };

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4, 
                        class T5> class F>
    struct quote5 {
      template <class T1, 
                class T2, 
                class T3,
                class T4,
                class T5> 
      struct apply : mpl::type_wrap< F<T1,T2,T3,T4,T5> >{};
    };

    template <class F>
    struct quote : quote0<F> {};

    template < template<class T1> class F,
               class T1>
    struct quote< F<T1> > : quote1<F> {};

    template < template<class T1, 
                        class T2> class F,
               class T1,
               class T2>
    struct quote< F<T1,T2> > : quote2<F> {};

    template < template<class T1, 
                        class T2, 
                        class T3> class F,
               class T1,
               class T2,
               class T3>
    struct quote< F<T1,T2,T3> > : quote3<F> {};

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4> class F,
               class T1,
               class T2,
               class T3,
               class T4>
    struct quote< F<T1,T2,T3,T4> > : quote4<F> {};

    template < template<class T1, 
                        class T2, 
                        class T3, 
                        class T4, 
                        class T5> class F,
               class T1,
               class T2,
               class T3,
               class T4,
               class T5>
    struct quote< F<T1,T2,T3,T4,T5> > : quote5<F> {};

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_QUOTE_HPP
