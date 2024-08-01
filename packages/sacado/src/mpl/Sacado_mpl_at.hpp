// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_AT_HPP
#define SACADO_MPL_AT_HPP

namespace Sacado {

  namespace mpl {

    template <class T, int Pos> struct at_impl {};

    template <class T, int Pos>
    struct at : 
      at_impl<typename T::tag,Pos>:: template apply<T> {};

  }

}

#endif // SACADO_MPL_AT_HPP
