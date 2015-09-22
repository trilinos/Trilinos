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

#ifndef __TSQR_TestSetup_hpp
#define __TSQR_TestSetup_hpp

#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_Random_GlobalMatrix.hpp>
#include <Tsqr_Matrix.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <vector>


namespace TSQR { 
  namespace Test {

    template<class Ordinal, class CommOrdinal>
    Ordinal 
    numLocalRows (const Ordinal nrowsGlobal, 
		  const CommOrdinal myRank,
		  const CommOrdinal nprocs)
    {
      const Ordinal nrowsLocal = nrowsGlobal / Ordinal(nprocs);
      const Ordinal remainder = nrowsGlobal - nrowsLocal * Ordinal(nprocs);
      if (myRank != nprocs - 1)
	return nrowsLocal;
      else
	return nrowsLocal + remainder;
    }

    /// \param generator [in/out] Proc 0 is the only MPI process that
    /// generates pseudorandom numbers.  This allows us to use a
    /// sequential PRNG.  Otherwise, we would have to use a parallel
    /// PRNG, in order to prevent correlations between numbers
    /// generated on different MPI processes.  On processes other than
    /// Proc 0, generator is not touched.
    template<class MatrixViewType, class Generator>
    void
    distributedTestProblem (Generator& generator,
			    MatrixViewType& A_local,
			    MessengerBase<typename MatrixViewType::ordinal_type>* const ordinalComm,
			    MessengerBase<typename MatrixViewType::scalar_type>* const scalarComm)
    {
      typedef typename MatrixViewType::ordinal_type ordinal_type;
      typedef typename MatrixViewType::scalar_type scalar_type;
      typedef typename Teuchos::ScalarTraits< scalar_type >::magnitudeType magnitude_type;

      const int myRank = scalarComm->rank();
      const ordinal_type ncols = A_local.ncols();

      if (myRank == 0)
	{
	  // Generate some singular values for the test problem.
	  std::vector< magnitude_type > singular_values (ncols);
	  singular_values[0] = 1.0;
	  for (ordinal_type k = 1; k < ncols; ++k)
	    singular_values[k] = singular_values[k-1] / double(2);

	  // Generate the test problem.  All MPI processes
	  // participate, but only Proc 0 generates the (pseudo)random
	  // numbers.
	  TSQR::Random::randomGlobalMatrix (&generator, A_local, 
					    &singular_values[0], ordinalComm,
					    scalarComm);
	}
      else
	{
	  // This helps C++ deduce the type; the values aren't read on
	  // this proc.
	  magnitude_type singular_values[1];

	  // All MPI processes participate in the distribution of the
	  // test matrix.
	  TSQR::Random::randomGlobalMatrix (&generator, A_local, 
					    &singular_values[0], ordinalComm,
					    scalarComm);
	}
    }
  } // namespace Test
} // namespace TSQR

#endif // __TSQR_TestSetup_hpp
