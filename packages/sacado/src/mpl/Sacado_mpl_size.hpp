// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_SIZE_HPP
#define SACADO_MPL_SIZE_HPP

namespace Sacado {

  namespace mpl {

    template <class T> struct size_impl {};

    template <class T>
    struct size : 
      size_impl<typename T::tag>:: template apply<T> {};

  }

}

#endif // SACADO_MPL_SIZE_HPP
