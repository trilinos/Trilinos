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
   \file   Amesos2_NewSolver_FunctionMap.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Wed Jul 27 12:53:10 MDT 2011

   \brief  Template for providing a mechanism to map function calls to the
           correct Solver function based on the scalar type of Matrices and
           MultiVectors
*/

#ifndef AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
#define AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
//#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_PardisoMKL_TypeMap.hpp"


namespace Amesos2 {

  namespace PMKL {
#   include "mkl_pardiso.h"
  }

  /** \internal
   * 
   * For Pardiso we bind to the library functions based on the local
   * ordinal type.  If the local ordinal type is bigger than int, then
   * we use pardiso_64 instead.  The void* arrays are interpreted by
   * the function based on the value of mtype and iparm(28) as
   * single/double and complex/real.
   */
  template <>
  struct FunctionMap<PardisoMKL,PMKL::_INTEGER_t>
  {
    static void pardiso( void* pt,
			 PMKL::_INTEGER_t* maxfct, PMKL::_INTEGER_t* mnum,
			 PMKL::_INTEGER_t* mtype , PMKL::_INTEGER_t* phase,
			 PMKL::_INTEGER_t* n     , void* a, PMKL::_INTEGER_t* ia,
			 PMKL::_INTEGER_t* ja    , PMKL::_INTEGER_t* perm,
			 PMKL::_INTEGER_t* nrhs  , PMKL::_INTEGER_t* iparm,
			 PMKL::_INTEGER_t* msglvl, void* b, void* x,
			 PMKL::_INTEGER_t* error)
    {
      PMKL::pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
		    perm, nrhs, iparm, msglvl, b, x, error);
    }
  };


  template <>
  struct FunctionMap<PardisoMKL,long long int>
  {
    static void pardiso( void* pt,
			 long long int* maxfct, long long int* mnum,
			 long long int* mtype , long long int* phase,
			 long long int* n     , void* a, long long int* ia,
			 long long int* ja    , long long int* perm,
			 long long int* nrhs  , long long int* iparm,
			 long long int* msglvl, void* b, void* x,
			 long long int* error)
    {
      PMKL::pardiso_64(pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
		       perm, nrhs, iparm, msglvl, b, x, error);
    }
  };

} // end namespace Amesos2

#endif  // AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
