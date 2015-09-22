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
   \file   Amesos2_FunctionMap.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Thu Jul 15 22:04:41 2010

   \brief  Declaration of Function mapping class for Amesos2.
*/


#ifndef AMESOS2_FUNCTIONMAP_HPP
#define AMESOS2_FUNCTIONMAP_HPP

namespace Amesos2 {

/**
 * \brief Passes functions to TPL functions based on type.
 *
 * Helper class which passes on function calls to the appropriate
 * Solver function based on the type of its scalar template argument.
 *
 * Some Solvers have solver and matrix builder functions defined based
 * on data type.  One function for complex, one for double precision
 * complex, another for \c float , and yet another for \c double.  To
 * work elegantly with the Amesos2::Solver interface we want to be
 * able to perform a single function call which is appropriate for the
 * scalar type of the Matrix and MultiVectors that we are working
 * with.  The \c FunctionMap class provides that capability.
 *
 * The class template is specialized for each Solver and each data
 * type that it supports, and errors are thrown for other data types.
 */
template <template <typename,typename> class ConcreteSolver, typename Scalar>
struct FunctionMap
{};


} // end namespace Amesos2

#endif  // AMESOS2_FUNCTIONMAP_HPP
