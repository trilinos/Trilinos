// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_IF_HPP
#define SACADO_MPL_IF_HPP

#include "Sacado_mpl_type_wrap.hpp"

namespace Sacado {

  namespace mpl {

    template <bool cond, class T1, class T2> struct mpl_if_c {};
    template <class T1, class T2> struct mpl_if_c<true,T1,T2> :
      type_wrap<T1> {};
    template <class T1, class T2> struct mpl_if_c<false,T1,T2> :
      type_wrap<T2> {};

    template <class C, class T1, class T2> struct mpl_if :
      mpl_if_c<C::value,T1,T2> {};

  }

}

#endif // SACADO_MPL_IF_HPP
