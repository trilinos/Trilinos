// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziTraceMinTypes.hpp
 *  Defines some enumerated types for TraceMin's parameters
*/

#ifndef ANASAZI_TRACEMIN_TYPES_HPP
#define ANASAZI_TRACEMIN_TYPES_HPP

namespace Anasazi {
namespace Experimental {

  enum WhenToShiftType
  {
    NEVER_SHIFT,
    SHIFT_WHEN_TRACE_LEVELS,
    SHIFT_WHEN_RESID_SMALL,
    ALWAYS_SHIFT
  };

  enum HowToShiftType
  {
    LARGEST_CONVERGED_SHIFT,
    ADJUSTED_RITZ_SHIFT,
    RITZ_VALUES_SHIFT,
    EXPERIMENTAL_SHIFT
  };

  enum SaddleSolType
  {
    PROJECTED_KRYLOV_SOLVER,
    SCHUR_COMPLEMENT_SOLVER,
    BD_PREC_MINRES,
    HSS_PREC_GMRES
  };

}}
#endif
