// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SCADO_MPL_VECTOR_SIZE_SPEC_HPP
#define SCADO_MPL_VECTOR_SIZE_SPEC_HPP

namespace Sacado {

  namespace mpl {

    template <typename... Args>
    struct vector_size {
      static const int sz = sizeof...(Args);
    };

  }

}

#endif // SCADO_MPL_VECTOR_SIZE_SPEC_HPP
