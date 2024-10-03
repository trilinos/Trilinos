// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_INTEGRAL_C_HPP
#define SACADO_MPL_INTEGRAL_C_HPP

namespace Sacado {

  namespace mpl {

    // Type wrapper for storing an integral value
    template <class T, T N>
    struct integral_c {
      static const T value = N;
      typedef integral_c<T,N> type;
      typedef T value_type;
      typedef integral_c<T,N+1> next;
      typedef integral_c<T,N-1> prior;
      operator T() const { return N; }
    };

  }

}

#endif
