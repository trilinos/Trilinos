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

#ifndef __Tsqr_CombineBenchmarker_hpp
#define __Tsqr_CombineBenchmarker_hpp

#include "Tsqr_CombineBenchmark.hpp"

#include <Tsqr_Random_NormalGenerator.hpp>
#include <Tsqr_Random_MatrixGenerator.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

#include <Tsqr_ApplyType.hpp>
#include <Tsqr_Matrix.hpp>
#include <Tsqr_ScalarTraits.hpp>
#include <Tsqr_Util.hpp>

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    template< class Ordinal, class Scalar, class CombineType, class TimerType >
    class CombineBenchmarker {
    public:
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;
      typedef TSQR::Random::NormalGenerator< Ordinal, Scalar > normgen_type;
      typedef TSQR::Random::NormalGenerator< Ordinal, magnitude_type > normgen_mag_type;
      typedef TSQR::Random::MatrixGenerator< Ordinal, Scalar, normgen_type > matgen_type;
      typedef Matrix< Ordinal, Scalar > matrix_type;
      typedef std::pair< magnitude_type, magnitude_type > results_type;

    private:
      normgen_type gen_;
      normgen_mag_type magGen_;
      int numTrials_;
      bool debug_;

    public:

      CombineBenchmarker (normgen_type& gen, 
			  normgen_mag_type& magGen, 
			  const int numTrials,
			  const bool debug = false) : 
	gen_ (gen), magGen_ (magGen), numTrials_ (numTrials), debug_ (debug)
      {
	TSQR::Test::verifyTimerConcept< TimerType >();
      }

      double 
      R1R2_benchmark (const Ordinal numCols)
      {
	using std::cerr;
	using std::endl;
	using std::invalid_argument;
	using std::ostringstream;
	using std::vector;

	if (numCols == 0)
	  throw invalid_argument ("ncols == 0 is not allowed");

	// Generate two different sets of singular values.
	// Randomly perturb them, but make sure all are positive.

	vector< magnitude_type > sigma_R1 (numCols);
	vector< magnitude_type > sigma_R2 (numCols);
	generateSingularValues (sigma_R1, numCols);
	generateSingularValues (sigma_R2, numCols);

	matrix_type R1 (numCols, numCols, Scalar(0));
	matrix_type R2 (numCols, numCols, Scalar(0));

	matgen_type matgen (gen_);
	matgen.fill_random_R (numCols, R1.get(), R1.lda(), &sigma_R1[0]);
	matgen.fill_random_R (numCols, R2.get(), R2.lda(), &sigma_R2[0]);

	// Space to put the original test problem, expressed as one
	// dense matrix rather than in two blocks.  This will be a
	// deep copy of the test problem, since the test problem
	// matrices will be overwritten by the factorizations.
	matrix_type A_R1R2 (Ordinal(2) * numCols, numCols, Scalar(0));

	// Copy [R1; R2] into A_R1R2.
	copy_matrix (numCols, numCols, &A_R1R2(0, 0), A_R1R2.lda(), 
		     R1.get(), R1.lda());
	copy_matrix (numCols, numCols, &A_R1R2(numCols, 0), A_R1R2.lda(), 
		     R2.get(), R2.lda());

	// Space to put the explicit Q factor
	matrix_type Q_R1R2 (Ordinal(2) * numCols, numCols, Scalar(0));

	// Fill the explicit Q factor matrix with the first numCols
	// columns of the identity matrix.
	for (Ordinal k = 0; k < numCols; ++k)
	  Q_R1R2(k, k) = Scalar(1);

	// tau factor array
	vector< Scalar > tau_R1R2 (numCols);

	// Workspace array for factorization and applying the Q factor.
	vector< Scalar > work (numCols);

	if (debug_)
	  cerr << endl << "----------------------------------------" << endl
	       << "TSQR::Combine test problem:" << endl
	       << "qr( [R1; R2] ), with R1 and R2 " << numCols 
	       << " by " << numCols << endl << endl;

	CombineType combiner; 

	// A few warmup runs just to avoid timing anomalies
	const int numWarmupRuns = 5;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_pair (numCols, R1.get(), R1.lda(), 
				  R2.get(), R2.lda(),
				  &tau_R1R2[0], &work[0]);
	    combiner.apply_pair (ApplyType("N"), numCols, numCols, 
				 R2.get(), R2.lda(), &tau_R1R2[0], 
				 &Q_R1R2(0, 0), Q_R1R2.lda(),
				 &Q_R1R2(numCols, 0), Q_R1R2.lda(),
				 &work[0]);
	  }

	// The actual benchmark timing runs
	TimerType timer ("Combine");
	timer.start ();
	for (int trialNum = 0; trialNum < numTrials_; ++trialNum)
	  {
	    combiner.factor_pair (numCols, R1.get(), R1.lda(), R2.get(), R2.lda(),
				  &tau_R1R2[0], &work[0]);
	    combiner.apply_pair (ApplyType("N"), numCols, numCols, 
				 R2.get(), R2.lda(), &tau_R1R2[0], 
				 &Q_R1R2(0, 0), Q_R1R2.lda(),
				 &Q_R1R2(numCols, 0), Q_R1R2.lda(),
				 &work[0]);
	  }
	const double timingResult = timer.stop();
	if (debug_)
	  cerr << "Results: " << numTrials_ << " consecutive runs took a total of " 
	       << timingResult << " seconds" << endl;
	return timingResult;
      }

      double 
      RA_benchmark (const Ordinal numRows, const Ordinal numCols)
      {
	using std::cerr;
	using std::endl;
	using std::invalid_argument;
	using std::ostringstream;
	using std::vector;

	if (numRows < numCols)
	  {
	    ostringstream os;
	    os << "# rows < # columns is not allowed.  You specified # rows = " 
	       << numRows << " and # columns = " << numCols << ".";
	    throw invalid_argument (os.str());
	  }
	else if (numCols == 0)
	  throw invalid_argument ("ncols == 0 is not allowed");

	// Generate two different sets of singular values.  Randomly
	// perturb them, but make sure all are positive.
	vector< magnitude_type > sigma_R3 (numCols);
	vector< magnitude_type > sigma_A (numCols);
	generateSingularValues (sigma_R3, numCols);
	generateSingularValues (sigma_A, numCols);

	// Generate the test problem [R3; A]
	matrix_type R3 (numCols, numCols, Scalar(0));
	matrix_type A (numRows, numCols, Scalar(0));
	matgen_type matgen (gen_);
	matgen.fill_random_R (numCols, R3.get(), R3.lda(), &sigma_R3[0]);
	matgen.fill_random_svd (numRows, numCols, A.get(), A.lda(), &sigma_A[0]);

	// Space to put the original test problem, expressed as one
	// dense matrix rather than in two blocks.  This will be a deep
	// copy of the test problem, since the test problem matrices
	// will be overwritten by the factorization.
	matrix_type A_R3A (numRows + numCols, numCols, Scalar(0));

	// Copy [R3; A] into A_R3A.
	copy_matrix (numCols, numCols, &A_R3A(0, 0), A_R3A.lda(), 
		     R3.get(), R3.lda());
	copy_matrix (numRows, numCols, &A_R3A(numCols, 0), A_R3A.lda(), 
		     A.get(), A.lda());

	// Space to put the explicit Q factor
	matrix_type Q_R3A (numRows + numCols, numCols, Scalar(0));

	// Fill the explicit Q factor matrix with the first numCols
	// columns of the identity matrix.
	for (Ordinal k = 0; k < numCols; ++k)
	  Q_R3A(k, k) = Scalar(1);

	// tau factor array
	vector< Scalar > tau_R3A (numCols);

	// Workspace array for factorization and applying the Q factor.
	vector< Scalar > work (numCols);

	if (debug_)
	  cerr << endl << "----------------------------------------" << endl
	       << "TSQR::Combine test problem:" << endl
	       << "qr( [R3; A] ), with R3 " << numCols << " by " << numCols 
	       << " and A " << numRows << " by " << numCols << endl << endl;

	CombineType combiner; 

	// A few warmup runs just to avoid timing anomalies
	const int numWarmupRuns = 5;
	for (int warmupRun = 0; warmupRun < numWarmupRuns; ++warmupRun)
	  {
	    combiner.factor_inner (numRows, numCols, R3.get(), R3.lda(),
				   A.get(), A.lda(), &tau_R3A[0], &work[0]);
	    combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau_R3A[0], 
				  &Q_R3A(0, 0), Q_R3A.lda(),
				  &Q_R3A(numCols, 0), Q_R3A.lda(), 
				  &work[0]);
	  }

	// The actual benchmark timing runs
	TimerType timer ("Combine");
	timer.start ();
	for (int trialNum = 0; trialNum < numTrials_; ++trialNum)
	  {
	    combiner.factor_inner (numRows, numCols, R3.get(), R3.lda(),
				   A.get(), A.lda(), &tau_R3A[0], &work[0]);
	    combiner.apply_inner (ApplyType("N"), numRows, numCols, numCols,
				  A.get(), A.lda(), &tau_R3A[0], 
				  &Q_R3A(0, 0), Q_R3A.lda(),
				  &Q_R3A(numCols, 0), Q_R3A.lda(), 
				  &work[0]);
	  }
	const double timingResult = timer.stop();
	if (debug_)
	  cerr << "Results: " << numTrials_ << " consecutive runs took a total of " 
	       << timingResult << " seconds" << endl;
	return timingResult;
      }

    private:
      void
      generateSingularValues (std::vector< magnitude_type >& sigma,
			      const Ordinal numValues)
      {
	const magnitude_type machEps = 
	  std::numeric_limits< magnitude_type >::epsilon();
	sigma.resize (numValues);
    
	// Relative amount by which to perturb each singular value.  The
	// perturbation will be multiplied by a normal(0,1) pseudorandom
	// number drawn from magGen.
	const magnitude_type perturbationFactor = magnitude_type(10) * machEps;

	sigma[0] = magnitude_type (1);
	for (Ordinal k = 1; k < numValues; ++k)
	  {
	    const magnitude_type perturbation = perturbationFactor * magGen_();
	    const magnitude_type beforePerturb = sigma[k-1] / magnitude_type(2);
	    const magnitude_type candidate = beforePerturb + perturbation;

	    // If adding the perturbation to beforePerturb would result
	    // in a nonpositive number, subtract instead.
	    if (candidate <= magnitude_type(0))
	      sigma[k] = beforePerturb - perturbation;
	    else
	      sigma[k] = candidate;
	  }
      }
    };

  } // namespace Test
} // namespace TSQR

#endif // __Tsqr_CombineBenchmarker_hpp
