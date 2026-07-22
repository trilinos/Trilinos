// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_SOLVER_MANAGER_MP_VECTOR_HPP
#define BELOS_SOLVER_MANAGER_MP_VECTOR_HPP

// Forward declaration
namespace Sacado {
  namespace MP {
    template <class S> class Vector;
  }
}

namespace Belos {
  namespace Details{

    // Forward declaration
    template<class S> class LapackSupportsScalar;

    // Declare MP::Vector scalar type supports LAPACK
    //
    // This isn't really true, but allows use of this scalar type in limited
    // circumstances in Belos that require LAPACK (e.g., PseudoBlockCG).
    template<class S>
    class LapackSupportsScalar< Sacado::MP::Vector<S> > {
    public:
      const static bool value = true;
    };
  }
}

#endif
