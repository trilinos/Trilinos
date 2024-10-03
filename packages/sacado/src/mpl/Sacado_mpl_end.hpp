// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_END_HPP
#define SACADO_MPL_END_HPP

namespace Sacado {

  namespace mpl {

    template <class T> struct end_impl {};

    template <class T>
    struct end : 
      end_impl<typename T::tag>:: template apply<T> {};

  }

}

#endif // SACADO_MPL_END_HPP
