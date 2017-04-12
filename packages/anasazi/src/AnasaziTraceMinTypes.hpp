// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
