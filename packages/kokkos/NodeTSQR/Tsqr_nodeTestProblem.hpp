//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
