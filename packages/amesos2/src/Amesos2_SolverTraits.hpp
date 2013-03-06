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
 * \file   Amesos2_SolverTraits.hpp
 * \author Eric Bavier <etbavie@muenster.srn.sandia.gov>
 * \date   Mon Jul 18 16:10:34 2011
 * 
 * \brief  Provides access to interesting solver traits
 * 
 * 
 */

#ifndef AMESOS2_SOLVERTRAITS_HPP
#define AMESOS2_SOLVERTRAITS_HPP

#include "Amesos2_Meta.hpp"

namespace Amesos2 {

  /** \internal
   *  
   * \brief Provides traits about solvers.
   *
   * \internal Each concrete solver interface should specialize this
   * struct for that particular solver.
   */
  template <template <class,class> class ConcreteSolver>
  struct solver_traits {
    typedef Meta::nil_t supported_scalars;
  };

  
  ////////////////////////////
  // Related meta-functions //
  ////////////////////////////
  
  /** \internal
   *
   * \brief Check whether a solver supports a scalar type
   * 
   * A meta-function for checking a solver's support for a scalar.
   * The special-case is a Meta::nil_t list, which we are going to
   * interpret as the solver supporting any and all types.
   * If a solver's supported_scalars type-list is nil, this
   * meta-function returns true always.  Otherwise, it checks whether
   * the given type is in the solver's supported_scalars list.
   */

    /* SR: We will not use external initialization for the static const types.
     * Combined with template meta programming this fails in Intel compilers
     * 11-13. Moving all the initializations inside the declarations.
     */
  template <template <class,class> class ConcreteSolver,
	    typename Scalar>
  struct solver_supports_scalar {
    static const bool value =
              Meta::if_then_else<Meta::is_same<typename solver_traits<ConcreteSolver>::supported_scalars, Meta::nil_t>::value,
		       Meta::true_type,
		       Meta::type_list_contains<
			 typename solver_traits<ConcreteSolver>::supported_scalars,
			 Scalar> >::type::value;
  };

} // end namespace Amesos2

#endif	// AMESOS2_SOLVERTRAITS_HPP
