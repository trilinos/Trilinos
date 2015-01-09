// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
    RITZ_VALUES_SHIFT
  };

  enum SaddleSolType
  {
    PROJECTED_KRYLOV_SOLVER,
    SCHUR_COMPLEMENT_SOLVER,
    BD_PREC_MINRES
  };

}}
#endif
