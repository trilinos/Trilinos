// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
//
// @HEADER

/**
   \file   Amesos2_KLU2_FunctionMap.hpp
   \author Siva Rajamanickam <srajama@sandia.gov>

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_KLU2_FUNCTIONMAP_HPP
#define AMESOS2_KLU2_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_KLU2_TypeMap.hpp"


/* External definitions of the KLU2 functions
 *
 */
 // TODO
namespace KLU2 {
#include "klu2_defaults.hpp"
#include "klu2_analyze.hpp"
#include "klu2_factor.hpp"
#include "klu2_solve.hpp"
#include "klu2_free_symbolic.hpp"
#include "klu2_free_numeric.hpp"
} // end namespace KLU2


namespace Amesos2 {

  /* ==================== Specializations ====================
   *
   * \cond KLU2_function_specializations
   */

  /**
   * \brief Pass function calls to KLU2 based on data type.
   *
   * Helper class which passes on function calls to the appropriate
   * KLU2 function based on the type of its scalar template argument.
   *
   * KLU2 has solver and matrix builder functions defined based on
   * data type.  One function for complex, one for double precision
   * complex, another for \c float , and yet another for \c double.  To
   * work elegantly with the Amesos2::KLU2 interface we want to be
   * able to perform a single function call which is appropriate for the
   * scalar type of the Matrix and MultiVectors that we are working
   * with.  The \c FunctionMap class provides that capability.
   *
   * The class template is specialized for each data type that KLU2
   * supports.  The Amesos2::create function assures that an
   * unspecialized FunctionMap will never be called by the solver
   * interface.
   *
   */
  // TODO : Do we need the specializations for KLU2 ??


  /* \endcond KLU2_function_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_KLU2_FUNCTIONMAP_HPP
