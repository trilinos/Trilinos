// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_BIND_HPP
#define SACADO_MPL_BIND_HPP

#include <type_traits>

#include "Sacado_mpl_placeholders.hpp"
#include "Sacado_mpl_apply_wrap.hpp"

namespace Sacado {

  namespace mpl {

    template <int k, class F, class T1, class T2, class T3, class T4, class T5>
    struct hk { typedef F type; };

    template <int k, int N, class T1, class T2, class T3, class T4, class T5> 
    struct hk<k,arg<N>,T1,T2,T3,T4,T5> : 
      apply_wrap<arg<N>,T1,T2,T3,T4,T5> {};

    template <int k, class T1, class T2, class T3, class T4, class T5> 
    struct hk<k,arg<-1>,T1,T2,T3,T4,T5> : 
      apply_wrap<arg<k>,T1,T2,T3,T4,T5> {};

    template <class F, class T1>
    struct bind1 {
      template <class U1=mpl::none, 
                class U2=mpl::none, 
                class U3=mpl::none, 
                class U4=mpl::none, 
                class U5=mpl::none> 
      struct apply : 
        apply_wrap1<F,
                    typename hk<std::is_same<T1,placeholders::_>::value,
                                T1,
                                U1,U2,U3,U4,U5>::type> {};
    };

    template <class F, class T1, class T2>
    struct bind2 {
      template <class U1=mpl::none, 
                class U2=mpl::none, 
                class U3=mpl::none, 
                class U4=mpl::none, 
                class U5=mpl::none> 
      struct apply : 
        apply_wrap2<F,
                    typename hk<std::is_same<T1,placeholders::_>::value,
                                T1,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value,
                                T2,
                                U1,U2,U3,U4,U5>::type> {};
    };

    template <class F, class T1, class T2, class T3>
    struct bind3 {
      template <class U1=mpl::none, 
                class U2=mpl::none, 
                class U3=mpl::none, 
                class U4=mpl::none, 
                class U5=mpl::none> 
      struct apply : 
        apply_wrap3<F,
                    typename hk<std::is_same<T1,placeholders::_>::value,
                                T1,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value,
                                T2,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value+
                                std::is_same<T3,placeholders::_>::value,
                                T3,
                                U1,U2,U3,U4,U5>::type> {};
    };

    template <class F, class T1, class T2, class T3, class T4>
    struct bind4 {
      template <class U1=mpl::none, 
                class U2=mpl::none, 
                class U3=mpl::none, 
                class U4=mpl::none, 
                class U5=mpl::none> 
      struct apply : 
        apply_wrap4<F,
                    typename hk<std::is_same<T1,placeholders::_>::value,
                                T1,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value,
                                T2,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value+
                                std::is_same<T3,placeholders::_>::value,
                                T3,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value+
                                std::is_same<T3,placeholders::_>::value+
                                std::is_same<T4,placeholders::_>::value,
                                T4,
                                U1,U2,U3,U4,U5>::type> {};
    };

    template <class F, class T1, class T2, class T3, class T4, class T5>
    struct bind5 {
      template <class U1=mpl::none, 
                class U2=mpl::none, 
                class U3=mpl::none, 
                class U4=mpl::none, 
                class U5=mpl::none> 
      struct apply : 
        apply_wrap5<F,
                    typename hk<std::is_same<T1,placeholders::_>::value,
                                T1,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value,
                                T2,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value+
                                std::is_same<T3,placeholders::_>::value,
                                T3,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value+
                                std::is_same<T3,placeholders::_>::value+
                                std::is_same<T4,placeholders::_>::value,
                                T4,
                                U1,U2,U3,U4,U5>::type,
                    typename hk<std::is_same<T1,placeholders::_>::value+
                                std::is_same<T2,placeholders::_>::value+
                                std::is_same<T3,placeholders::_>::value+
                                std::is_same<T4,placeholders::_>::value+
                                std::is_same<T5,placeholders::_>::value,
                                T5,
                                U1,U2,U3,U4,U5>::type> {};
    };

  } // namespace mpl

} // namespace Sacado

#endif // SACADO_MPL_BIND_HPP
