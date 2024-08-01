// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_IS_PLACEHOLDER_HPP
#define SACADO_MPL_IS_PLACEHOLDER_HPP

#include "Sacado_mpl_is_placeholder.hpp"
#include "Sacado_mpl_none.hpp"

namespace Sacado {

  namespace mpl {

    template <class F> 
    struct is_placeholder { 
      static const bool value = false; 
    };
    template <int N> 
    struct is_placeholder< arg<N> > {
      static const bool value = true; 
    };
    template <template <class T1> class F,
              class T1> 
    struct is_placeholder< F<T1> > {
      static const bool value = is_placeholder<T1>::value;
    };
    template <template <class T1, class T2> class F,
              class T1,
              class T2> 
    struct is_placeholder< F<T1,T2> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value;
    };
    template <template <class T1, class T2, class T3> class F,
              class T1,
              class T2,
              class T3> 
    struct is_placeholder< F<T1,T2,T3> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value || 
        is_placeholder<T3>::value;
    };
    template <template <class T1, class T2, class T3, class T4> class F,
              class T1,
              class T2,
              class T3,
              class T4> 
    struct is_placeholder< F<T1,T2,T3,T4> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value || 
        is_placeholder<T3>::value || 
        is_placeholder<T4>::value;
    };
    template <template <class T1, class T2, class T3, class T4, class T5> class F,
              class T1,
              class T2,
              class T3,
              class T4,
              class T5> 
    struct is_placeholder< F<T1,T2,T3,T4,T5> > {
      static const bool value = 
        is_placeholder<T1>::value || 
        is_placeholder<T2>::value || 
        is_placeholder<T3>::value || 
        is_placeholder<T4>::value || 
        is_placeholder<T5>::value;
    };

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_IS_PLACEHOLDER_HPP
