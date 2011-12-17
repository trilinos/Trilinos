//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __TSQR_Test_nodeTestProblem_hpp
#define __TSQR_Test_nodeTestProblem_hpp

#include <Tsqr_Random_MatrixGenerator.hpp>
#include <Tsqr_ScalarTraits.hpp>

#include <algorithm>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    /// Fill in the nrows by ncols matrix A (with leading dimension
    /// lda) with a test problem for single-node TSQR.
    template< class Ordinal, class Scalar, class Generator >
    void
    nodeTestProblem (Generator& generator,
		     const Ordinal nrows, 
		     const Ordinal ncols,
		     Scalar A[],
		     const Ordinal lda,
		     const bool numerically_interesting)
    {
      typedef TSQR::Random::MatrixGenerator< Ordinal, Scalar, Generator > matgen_type; 
      matgen_type matGen (generator);
      
      if (numerically_interesting)
	{
	  typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

	  std::vector< magnitude_type > singular_values (std::min(nrows, ncols));
	  singular_values[0] = magnitude_type (1);
	  for (Ordinal k = 1; k < std::min(nrows, ncols); ++k)
	    singular_values[k] = singular_values[k-1] / magnitude_type(2);

	  matGen.fill_random_svd (nrows, ncols, A, lda, &singular_values[0]);
	}
      else
	matGen.fill_random (nrows, ncols, A, lda);
    }

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_nodeTestProblem_hpp
