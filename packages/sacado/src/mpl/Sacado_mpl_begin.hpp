// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_MPL_BEGIN_HPP
#define SACADO_MPL_BEGIN_HPP

namespace Sacado {

  namespace mpl {

    template <class T> struct begin_impl {};

    template <class T>
    struct begin : 
      begin_impl<typename T::tag>:: template apply<T> {};

  }

}

#endif // SACADO_MPL_BEGIN_HPP
