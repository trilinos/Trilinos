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
