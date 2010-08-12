// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef __TSQR_Test_ParTest_hpp
#define __TSQR_Test_ParTest_hpp

#include <Tsqr_generateStack.hpp>
#include <Tsqr_MessengerBase.hpp>
#include <Tsqr_DistTsqr.hpp>
#include <Tsqr_GlobalVerify.hpp>
#include <Tsqr_printGlobalMatrix.hpp>
#include <Tsqr_ScalarTraits.hpp>

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    /// Verify DistParallelTSQR functions
    ///
    /// \param ncols [in] Number of columns in the matrix to test.
    ///   (Number of rows := (# MPI processors) * ncols
    /// \param ordinalComm [in/out] Communicator object over which to test.
    /// \param scalarComm [in/out] Communicator object over which to test.
    /// \param b_debug [in] Whether to show verbose debug output to stderr.
    template< class Ordinal, class Scalar, class Generator >
    void
    verifyParTsqr (Generator& generator,
		   const Ordinal ncols, 
		   MessengerBase< Ordinal >* const ordinalComm,
		   MessengerBase< Scalar >* const scalarComm,
		   const bool b_debug = false)
    {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      const bool b_extra_debug = false;
      const int my_rank = scalarComm->rank();
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "par_tsqr_verify:" << endl;
	  scalarComm->barrier();
	}
      const Ordinal nrows_local = ncols;

      // A_local: Space for the matrix A to factor -- local to each
      //   processor.
      // A_global: Global matrix (only nonempty on Proc 0)
      Matrix< Ordinal, Scalar > A_local, A_global;
      // This modifies A_local on all procs, and A_global as well on
      // Proc 0.
      par_tsqr_test_problem (generator, A_local, A_global, ncols, scalarComm);

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Generated test problem." << endl;
	  scalarComm->barrier();
	}

      // Copy the test problem input into R, since the factorization will
      // overwrite it place with the final R factor.
      Matrix< Ordinal, Scalar > R (ncols, ncols);
      R.copy (A_local);

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished copying test problem input into (local) R." << endl;
	}

      // Make sure that the test problem (the matrix to factor) was
      // distributed correctly.
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- R stack test matrix:" << endl;
	  scalarComm->barrier();
	  printGlobalMatrix (cerr, A_local, scalarComm, ordinalComm);
	}

      // Set up TSQR.
      DistTsqr< Ordinal, Scalar > par (scalarComm);
      if (b_extra_debug && b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- All MPI process(es) have initialized their "
	      "DistTsqr object." << endl << endl;
	}

      // Factor the matrix A (copied into R, which will be overwritten on output)
      typedef typename DistTsqr< Ordinal, Scalar >::FactorOutput factor_output_type;
      factor_output_type factor_output = par.factor (R.view());

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished DistTsqr::factor" << endl;
	}

      // Prepare space in which to construct the explicit Q factor (local
      // component on this processor)
      Matrix< Ordinal, Scalar > Q_local (nrows_local, ncols);

      // Compute the explicit Q factor
      par.explicit_Q (ncols, Q_local.get(), Q_local.lda(), factor_output);

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished DistTsqr::explicit_Q" << endl;
	}

      // Verify the factorization
      std::pair< magnitude_type, magnitude_type > result = 
	global_verify (nrows_local, ncols, A_local.get(), A_local.lda(),
		       Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
		       scalarComm);
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished global_verify" << endl;
	}

      scalarComm->barrier();
      if (my_rank == 0)
	{
	  cout << "Relative residual $\\|A - Q*R\\|_2 / \\|A\\|_2$ = " 
	       << result.first << endl;
	  cout << "Relative orthogonality $\\|I - Q^T*Q\\|_2$ = " 
	       << result.second << endl;
	}
    }

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_ParTest_hpp
