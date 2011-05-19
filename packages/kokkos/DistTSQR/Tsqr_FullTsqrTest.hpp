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

#ifndef __TSQR_Test_FullTsqrTest_hpp
#define __TSQR_Test_FullTsqrTest_hpp

#include <Tsqr.hpp>
#include <Tsqr_Random_NormalGenerator.hpp>
#include <Tsqr_Random_GlobalMatrix.hpp>
#include <Tsqr_TestSetup.hpp>
//#include <TsqrFactory_SequentialTsqr.hpp>
#include <Tsqr_GlobalVerify.hpp>
#include <Tsqr_TeuchosMessenger.hpp>
#include "Tsqr_TestUtils.hpp"
#include <Teuchos_ScalarTraits.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

namespace TSQR {
  namespace Test {

    /// \class TsqrInaccurate
    /// \brief Signals that a TSQR test failed due to insufficient accuracy.
    ///
    class TsqrInaccurate : public std::exception {
    public:
      //! Constructor
      TsqrInaccurate (const std::string& msg) : msg_ (msg) {}

      //! The error message
      const char* what() const throw() { return msg_.c_str(); }

      //! Destructor (declared virtual for memory safety of subclasses).
      virtual ~TsqrInaccurate() throw() {}

    private:
      std::string msg_;
    };

    /// \class FullTsqrVerifier
    /// \brief Test (correctness and) accuracy of Tsqr for one Scalar type.
    /// \author Mark Hoemmen
    ///
    /// This class is meant to be used only by \c
    /// FullTsqrVerifierCaller.  It performs one accuracy test of \c
    /// Tsqr for the given Scalar type (that is, the type of the
    /// matrix entries).  An accuracy test is also a correctness test.
    /// This test computes accuracy bounds for both orthogonality and
    /// forward errors, and if those bounds are exceeded and the
    /// failIfInaccurate option is enabled, the test will throw a \c
    /// TsqrInaccurate exception.
    ///
    /// The test takes a \c Teuchos::ParameterList input.  For a
    /// ParameterList with all parameters, their default values, and
    /// documentation, see the relevant class method in \c
    /// FullTsqrVerifierCaller.
    ///
    /// This class currently only tests the version of Tsqr that is
    /// the composition of NodeTsqrType=SequentialTsqr and
    /// DistTsqrType=DistTsqr.  This should suffice to test
    /// correctness, as long as the other NodeTsqrType possibilities
    /// (such as TbbTsqr) are tested separately.
    ///
    template<class Scalar>
    class FullTsqrVerifier {
    public:
      typedef Scalar scalar_type;
      typedef int ordinal_type;
      typedef SequentialTsqr<ordinal_type, scalar_type> node_tsqr_type;
      typedef DistTsqr<ordinal_type, scalar_type> dist_tsqr_type;
      typedef Tsqr<ordinal_type, scalar_type, node_tsqr_type> tsqr_type;
      typedef Kokkos::SerialNode node_type;

    private:

      //! Instantiate and return a (full) Tsqr instance.
      static Teuchos::RCP<tsqr_type>
      getTsqr (const Teuchos::RCP<Teuchos::ParameterList>& testParams, 
	       const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
      {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::rcp_implicit_cast;
	using Teuchos::RCP;
	using Teuchos::rcp;

	const size_t cacheSizeHint = testParams->get<size_t> ("cacheSizeHint");
	//const int numCores = testParams->get<int> ("numCores");

	//RCP<ParameterList> tsqrParams = parameterList ("Intranode TSQR");
	//tsqrParams->set ("cacheSizeHint", cacheSizeHint);
	//tsqrParams->set ("numCores", numCores);

	RCP<node_tsqr_type> seqTsqr = rcp (new node_tsqr_type (cacheSizeHint));

	RCP<TeuchosMessenger<scalar_type> > scalarMess =
	  rcp (new TeuchosMessenger<scalar_type> (comm));
	RCP<MessengerBase<scalar_type> > scalarMessBase = 
	  rcp_implicit_cast<MessengerBase<scalar_type> > (scalarMess);
	RCP<dist_tsqr_type> distTsqr = 
	  rcp (new dist_tsqr_type (scalarMessBase));

	return rcp (new tsqr_type (seqTsqr, distTsqr));
      }

    public:

      /// \brief Run the test for the Scalar type.
      ///
      /// \param comm [in] Communicator over which to run the test.
      /// \param node [in] Kokkos Node instance.
      /// \param testParams [in/out] Parameters for the test.  May
      ///   be modified by each test in turn.
      /// \param randomSeed [in/out] On input: the random seed for
      ///   LAPACK's pseudorandom number generator.  On output: the
      ///   updated random seed.
      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
	   const Teuchos::RCP<const node_type>& node,
	   const Teuchos::RCP<Teuchos::ParameterList>& testParams,
	   std::vector<int>& randomSeed)
      {
	using std::cerr;
	using std::cout;
	using std::endl;
	using Teuchos::arcp;
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::rcp_const_cast;
	using Teuchos::rcp_implicit_cast;
	typedef Matrix<ordinal_type, scalar_type> matrix_type;
	typedef MatView<ordinal_type, scalar_type> view_type;
	typedef typename tsqr_type::FactorOutput factor_output_type;

	const int myRank = Teuchos::rank (*comm);
	const int numProcs = Teuchos::size (*comm);

	// Construct TSQR implementation instance.
	RCP<tsqr_type> tsqr = getTsqr (testParams, comm);

	// Fetch test parameters from the input parameter list.
	const ordinal_type numRowsLocal = testParams->get<ordinal_type> ("numRowsLocal");
	const ordinal_type numCols = testParams->get<ordinal_type> ("numCols");	
	const int numCores = testParams->get<int> ("numCores");
	const bool contiguousCacheBlocks = testParams->get<bool> ("contiguousCacheBlocks");
	const bool testFactorExplicit = testParams->get<bool> ("testFactorExplicit");
	const bool testRankRevealing = testParams->get<bool> ("testRankRevealing");
	const bool debug = testParams->get<bool> ("debug");

	// Space for each node's local part of the test problem.
	// A_local, A_copy, and Q_local are distributed matrices, and
	// R is replicated on all processes sharing the communicator.
	matrix_type A_local (numRowsLocal, numCols);
	matrix_type A_copy (numRowsLocal, numCols);
	matrix_type Q_local (numRowsLocal, numCols);
	matrix_type R (numCols, numCols);

	// Start out by filling the test problem with zeros.
	typedef Teuchos::ScalarTraits<scalar_type> STS;
	A_local.fill (STS::zero());
	A_copy.fill (STS::zero());
	Q_local.fill (STS::zero());
	R.fill (STS::zero());

	// Create some reasonable singular values for the test problem:
	// 1, 1/2, 1/4, 1/8, ...
	typedef typename STS::magnitudeType magnitude_type;
	std::vector<magnitude_type> singularValues (numCols);
	typedef Teuchos::ScalarTraits<magnitude_type> STM;
	{
	  const magnitude_type scalingFactor = STM::one() + STM::one();
	  magnitude_type curVal = STM::one();
	  typedef typename std::vector<magnitude_type>::iterator iter_type;
	  for (iter_type it = singularValues.begin(); 
	       it != singularValues.end(); ++it)
	    {
	      *it = curVal;
	      curVal = curVal / scalingFactor;
	    }
	}

	// Construct a normal(0,1) pseudorandom number generator with
	// the given random seed.
	using TSQR::Random::NormalGenerator;
	typedef NormalGenerator<ordinal_type, scalar_type> generator_type;
	generator_type gen (randomSeed);

	// We need a Messenger for Ordinal-type data, so that we can
	// build a global random test matrix.
	RCP<MessengerBase<ordinal_type> > ordinalMessenger = 
	  rcp_implicit_cast<MessengerBase<ordinal_type> > (rcp (new TeuchosMessenger<ordinal_type> (comm)));

	// We also need a Messenger for Scalar-type data.  The TSQR
	// implementation already constructed one, but it's OK to
	// construct another one; TeuchosMessenger is just a thin
	// wrapper over the Teuchos::Comm object.
	RCP<MessengerBase<scalar_type> > scalarMessenger =
	  rcp_implicit_cast<MessengerBase<scalar_type> > (rcp (new TeuchosMessenger<scalar_type> (comm)));

	{
	  // Generate a global distributed matrix (whose part local to
	  // this node is in A_local) with the given singular values.
	  // This part has O(P) communication for P MPI processes.
	  using TSQR::Random::randomGlobalMatrix;
	  // Help the C++ compiler with type inference.
	  view_type A_local_view (A_local.nrows(), A_local.ncols(), A_local.get(), A_local.lda());
	  const magnitude_type* const singVals = (numCols == 0) ? NULL : &singularValues[0];
	  randomGlobalMatrix<view_type, generator_type> (&gen, A_local_view, singVals,
							 ordinalMessenger.getRawPtr(),
							 scalarMessenger.getRawPtr());
	}
	// Save the pseudorandom number generator's seed for any later
	// tests.  The generator keeps its own copy of the seed and
	// updates it internally, so we have to ask for its copy.
	gen.getSeed (randomSeed);

	// If specified in the test parameters, rearrange cache blocks
	// in the copy.  Otherwise, just copy the test problem into
	// A_copy.  The factorization overwrites the input matrix, so
	// we have to make a copy in order to validate the final
	// result.
	if (contiguousCacheBlocks)
	  {
	    tsqr->cache_block (numRowsLocal, numCols, A_copy.get(), 
			       A_local.get(), A_local.lda());
	    if (debug)
	      {
		Teuchos::barrier (*comm);
		if (myRank == 0)
		  cerr << "-- Finished Tsqr::cache_block" << endl;
	      }
	  }
	else
	  A_copy.copy (A_local);

	// "factorExplicit" is an alternate, hopefully faster way of
	// factoring the matrix, when only the explicit Q factor is
	// wanted.
	typedef Kokkos::MultiVector<scalar_type, node_type> KMV;
	if (testFactorExplicit)
	  {
	    // Kokkos::MultiVector wants a non-const Node, for some reason.
	    KMV A_copy_view (rcp_const_cast<node_type> (node));
	    A_copy_view.initializeValues (static_cast<size_t> (A_copy.nrows()), 
					  static_cast<size_t> (A_copy.ncols()), 
					  arcp (A_copy.get(), 0, A_copy.nrows()*A_copy.ncols(), false), // non-owning ArrayRCP
					  static_cast<size_t> (A_copy.lda()));
	    // Kokkos::MultiVector wants a non-const Node, for some reason.
	    KMV Q_view (rcp_const_cast<node_type> (node));
	    Q_view.initializeValues (static_cast<size_t> (Q_local.nrows()), 
				     static_cast<size_t> (Q_local.ncols()),
				     arcp (Q_local.get(), 0, Q_local.nrows()*Q_local.ncols(), false), // non-owning ArrayRCP
				     static_cast<size_t> (Q_local.lda()));
	    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> 
	      R_view (Teuchos::View, R.get(), R.lda(), R.nrows(), R.ncols());

	    tsqr->factorExplicit (A_copy_view, Q_view, R_view,
				  contiguousCacheBlocks);
	    if (debug)
	      {
		Teuchos::barrier (*comm);
		if (myRank == 0)
		  cerr << "-- Finished Tsqr::factorExplicit" << endl;
	      }
	  }
	else
	  {
	    // Factor the (copy of the) matrix.
	    factor_output_type factorOutput = 
	      tsqr->factor (numRowsLocal, numCols, A_copy.get(), A_copy.lda(), 
			    R.get(), R.lda(), contiguousCacheBlocks);
	    if (debug)
	      {
		Teuchos::barrier (*comm);
		if (myRank == 0)
		  cerr << "-- Finished Tsqr::factor" << endl;
	      }
	    // Compute the explicit Q factor in Q_local.
	    tsqr->explicit_Q (numRowsLocal, numCols, A_copy.get(), A_copy.lda(), 
			      factorOutput, numCols, Q_local.get(), Q_local.lda(),
			      contiguousCacheBlocks);
	    if (debug)
	      {
		Teuchos::barrier (*comm);
		if (myRank == 0)
		  cerr << "-- Finished Tsqr::explicit_Q" << endl;
	      }
	  }

	// Optionally, test rank-revealing capability.  We do this
	// before un-cache-blocking the explicit Q factor, since
	// revealRank can work with contiguous cache blocks, and
	// modifies the Q factor if the matrix doesn't have full
	// column rank.
	if (testRankRevealing)
	  {
	    // Kokkos::MultiVector wants a non-const Node, for some reason.
	    KMV Q_view (rcp_const_cast<node_type> (node));
	    Q_view.initializeValues (static_cast<size_t> (Q_local.nrows()), 
				     static_cast<size_t> (Q_local.ncols()),
				     arcp (Q_local.get(), 0, Q_local.nrows()*Q_local.ncols(), false), // non-owning ArrayRCP
				     static_cast<size_t> (Q_local.lda()));
	    Teuchos::SerialDenseMatrix<ordinal_type, scalar_type> 
	      R_view (Teuchos::View, R.get(), R.lda(), R.nrows(), R.ncols());
	    // If 2^{# columns} > machine precision, then our choice
	    // of singular values will make the smallest singular
	    // value < machine precision.  In that case, the SVD can't
	    // promise it will distinguish between tiny and zero.  If
	    // the number of columns is less than that, we can use a
	    // tolerance of zero to test the purported rank with the
	    // actual numerical rank.
	    const magnitude_type tol = STM::zero();
	    const ordinal_type rank = 
	      tsqr->revealRank (Q_view, R_view, tol, contiguousCacheBlocks);

	    magnitude_type two_to_the_numCols = STM::one();
	    for (int k = 0; k < numCols; ++k)
	      {
		const magnitude_type two = STM::one() + STM::one();
		two_to_the_numCols *= two;
	      }
	    // Throw in a factor of 10, just for more tolerance of
	    // rounding error (so the test only fails if something is
	    // really broken).
	    if (two_to_the_numCols > magnitude_type(10) * STM::eps())
	      {
		TEST_FOR_EXCEPTION(rank != numCols, std::logic_error,
				   "The matrix of " << numCols << " columns "
				   "should have full numerical rank, but Tsqr "
				   "reports that it has rank " << rank << ".  "
				   "Please report this bug to the Kokkos "
				   "developers.");
		if (debug)
		  {
		    Teuchos::barrier (*comm);
		    if (myRank == 0)
		      cerr << "-- Tested rank-revealing capability" << endl;
		  }
	      }
	    else
	      {
		if (debug)
		  {
		    Teuchos::barrier (*comm);
		    if (myRank == 0)
		      cerr << "-- Not testing rank-revealing capability; too many columns" << endl;
		  }
	      }
	  }
	// "Un"-cache-block the output, if contiguous cache blocks
	// were used.  This is only necessary because global_verify()
	// doesn't currently support contiguous cache blocks.
	if (contiguousCacheBlocks)
	  {
	    // We can use A_copy as scratch space for
	    // un-cache-blocking Q_local, since we're done using
	    // A_copy for other things.
	    tsqr->un_cache_block (numRowsLocal, numCols, A_copy.get(), 
				  A_copy.lda(), Q_local.get());
	    // Overwrite Q_local with the un-cache-blocked Q factor.
	    Q_local.copy (A_copy);
	    if (debug)
	      {
		Teuchos::barrier (*comm);
		if (myRank == 0)
		  cerr << "-- Finished Tsqr::un_cache_block" << endl;
	      }
	  }

	// Test accuracy of the factorization.
	const std::vector<magnitude_type> results = 
	  global_verify (numRowsLocal, numCols, A_local.get(), A_local.lda(),
			 Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
			 scalarMessenger.getRawPtr());
	if (debug)
	  {
	    Teuchos::barrier (*comm);
	    if (myRank == 0)
	      cerr << "-- Finished global_verify" << endl;
	  }

	// Print the results on Proc 0.
	if (myRank == 0)
	  {
	    if (testParams->get<bool> ("printFieldNames"))
	      {
		cout << "%"
		     << "method"
		     << ",scalarType"
		     << ",numRowsLocal"
		     << ",numCols"
		     << ",numProcs"
		     << ",numCores"
		     << ",cacheSizeHint"
		     << ",contiguousCacheBlocks"
		     << ",absFrobResid"
		     << ",absFrobOrthog"
		     << ",frobA" << endl;
		// We don't need to print field names again for the other tests, so set the test parameters accordingly.
		testParams->set ("printFieldNames", false);
	      }
	    if (testParams->get<bool> ("printResults"))
	      {
		cout << "Tsqr"
		     << "," << Teuchos::TypeNameTraits<scalar_type>::name()
		     << "," << numRowsLocal
		     << "," << numCols
		     << "," << numProcs
		     << "," << numCores
		     << "," << tsqr->cache_size_hint()
		     << "," << contiguousCacheBlocks 
		     << "," << results[0] 
		     << "," << results[1]
		     << "," << results[2]
		     << endl;
	      }
	  } // if (myRank == 0)

	// If requested, check accuracy and fail if results are not
	// sufficiently accurate.
	if (testParams->get<bool> ("failIfInaccurate"))
	  {
	    // Avoid overflow of the local Ordinal type, by casting
	    // first to a floating-point type.
	    const magnitude_type dimsProd = magnitude_type(numRowsLocal) * 
	      magnitude_type(numProcs) * magnitude_type(numCols*numCols);

	    // Relative residual error is ||A-Q*R|| / ||A||, or just
	    // ||A-Q*R|| if ||A|| == 0.  (The result had better be
	    // zero in the latter case.)  A reasonable error bound
	    // should incorporate the dimensions of the matrix, since
	    // this indicates the amount of rounding error.  Square
	    // root of the matrix dimensions is an old heuristic from
	    // Wilkinson or perhaps even an earlier source.  We
	    // include a factor of 10 so that the test won't fail
	    // unless there is a really good reason.
	    const magnitude_type relResidBound = 
	      magnitude_type(10) * STM::squareroot(dimsProd) * STM::eps();

	    // Orthogonality of the matrix should not depend on the
	    // matrix dimensions, if we measure in the 2-norm.
	    // However, we are measuring in the Frobenius norm, so
	    // it's appropriate to multiply eps by the number of
	    // entries in the matrix for which we compute the
	    // Frobenius norm.  We include a factor of 10 for the same
	    // reason as mentioned above.
	    const magnitude_type orthoBound = 
	      magnitude_type(10*numCols*numCols) * STM::eps();

	    // Avoid division by zero.
	    const magnitude_type relResidError = 
	      results[0] / (results[2] == STM::zero() ? STM::one() : results[2]);
	    TEST_FOR_EXCEPTION(relResidError > relResidBound, TsqrInaccurate,
			       "Full Tsqr (SequentialTsqr + DistTsqr) has an "
			       "inaccurate relative residual ||A - QR||_F"
			       << (results[2] == STM::zero() ? " / ||A||_F" : "")
			       << " = " << relResidError << ", which is greater"
			       " than the bound " << relResidBound << " by a "
			       "factor of " << relResidError / relResidBound
			       << ".");
	    const magnitude_type orthoError = results[1];
	    TEST_FOR_EXCEPTION(orthoError > orthoBound, TsqrInaccurate,
			       "Full Tsqr (SequentialTsqr + DistTsqr) has an "
			       "inaccurate orthogonality measure ||I - Q^* Q||"
			       "_F" << results[1] << " = " << orthoError 
			       << ", which is greater than the bound " 
			       << orthoBound << " by a factor of " 
			       << orthoError / orthoBound << ".");
	  } // if (the tests should fail on inaccuracy)
      }
    };

    /// \class FullTsqrVerifierCallerImpl
    /// \brief This class implements a "function template specialization."
    /// \author Mark Hoemmen
    ///
    /// We want to make FullTsqrVerifierCaller::run() a template
    /// function, with a partial specialization for Cons<CarType,
    /// CdrType> and a full specialization for NullType.  However,
    /// function templates can't have partial specializations, at
    /// least not in the version of the C++ standard currently
    /// supported by Trilinos.  Thus, I've taken the advice of Herb
    /// Sutter (C/C++ Users Journal, 19(7), July 2001), which can be
    /// read online here:
    ///
    /// http://www.gotw.ca/publications/mill17.htm
    ///
    /// Namely, I've implemented the function template via a class
    /// template.  This class is an implementation detail and not
    /// meant to be used anywhere else other than in
    /// FullTsqrVerifierCaller::run().
    template<class TypeListType>
    class FullTsqrVerifierCallerImpl {
    public:
      typedef Kokkos::SerialNode node_type; 

      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
	   const Teuchos::RCP<const node_type>& node,
	   const Teuchos::RCP<Teuchos::ParameterList>& testParams, 
	   std::vector<int>& randomSeed);
    };

    //
    // Partial specialization for Cons<CarType, CdrType>.
    //
    template<class CarType, class CdrType>
    class FullTsqrVerifierCallerImpl<TSQR::Test::Cons<CarType, CdrType> > {
    public:
      typedef Kokkos::SerialNode node_type; 

      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
	   const Teuchos::RCP<const node_type>& node,
	   const Teuchos::RCP<Teuchos::ParameterList>& testParams, 
	   std::vector<int>& randomSeed)
      {
	typedef CarType car_type;
	typedef CdrType cdr_type;
	FullTsqrVerifier<car_type>::run (comm, node, testParams, randomSeed);
	FullTsqrVerifierCallerImpl<cdr_type>::run (comm, node, testParams, randomSeed);
      }
    };

    //
    // Full specialization for NullCons.
    //
    template<>
    class FullTsqrVerifierCallerImpl<TSQR::Test::NullCons> {
    public:
      typedef Kokkos::SerialNode node_type; 

      static void
      run (const Teuchos::RCP<const Teuchos::Comm<int> >&,
	   const Teuchos::RCP<const node_type>&,
	   const Teuchos::RCP<Teuchos::ParameterList>&, 
	   std::vector<int>&)
      {
	// We're at the end of the type list, so do nothing.
      }
    };

    /// \class FullTsqrVerifierCaller
    /// \brief Invokes FullTsqrVerifier::run() over all Scalar types in a type list.
    /// \author Mark Hoemmen
    ///
    /// Use this class to test the full TSQR implementation in Tsqr.
    /// It will test Tsqr over a list of Scalar types that you define,
    /// using \c Cons and \c NullCons.
    class FullTsqrVerifierCaller {
    public:
      /// \typedef node_type
      /// \brief The Kokkos Node type to use.
      typedef Kokkos::SerialNode node_type;

      /// \typedef ordinal_type
      /// \brief The (local) Ordinal type to use for TSQR.
      ///
      /// This must be a type for which Teuchos::BLAS<ordinal_type,
      /// Scalar> and Teuchos::LAPACK<ordinal_type, Scalar> each have
      /// an instantiation.  That means a signed integer type.  LAPACK
      /// and the BLAS can be built with signed 64-bit integers
      /// (int64_t), but usually they are only built with signed
      /// 32-bit integers (int).
      typedef int ordinal_type;

      /// \brief Return a valid parameter list for verifying Tsqr.
      ///
      /// Call this once to get a valid parameter list with all the
      /// defaults filled in.  This list is valid for all the Scalar
      /// types which TsqrVerifierCaller::run tests.
      Teuchos::RCP<const Teuchos::ParameterList>
      getValidParameterList () const
      {
	using Teuchos::ParameterList;
	using Teuchos::parameterList;
	using Teuchos::RCP;

	RCP<ParameterList> plist = parameterList ("FullTsqrVerifier");

	const size_t cacheSizeHint = 0;
	const int numCores = 1;
	const ordinal_type numRowsLocal = 100;
	const ordinal_type numCols = 10;
	const bool contiguousCacheBlocks = false;
	const bool testFactorExplicit = true;
	const bool testRankRevealing = true;
	const bool printFieldNames = true;
	const bool printResults = true;
	const bool failIfInaccurate = true;
	const bool debug = false;

	// Parameters for configuring Tsqr itself.
	plist->set ("cacheSizeHint", cacheSizeHint, 
		    "Cache size hint in bytes.  "
		    "Zero means TSQR picks a reasonable default.");
	plist->set ("numCores", numCores,
		    "Number of partition(s) to use for TbbTsqr (if "
		    "applicable).  Must be a positive integer.");

	// Parameters for testing Tsqr.
	plist->set ("numRowsLocal", numRowsLocal, 
		    "Number of rows per (MPI) process in the test matrix.  "
		    "Must be >= the number of columns.");
	plist->set ("numCols", numCols, 
		    "Number of columns in the test matrix.");
	plist->set ("contiguousCacheBlocks", contiguousCacheBlocks, 
		    "Whether to test the factorization with contiguously "
		    "stored cache blocks.");
	plist->set ("testFactorExplicit", testFactorExplicit, 
		    "Whether to test TSQR's factorExplicit() (a hopefully "
		    "faster path than calling factor() and explicit_Q() in "
		    "sequence).");
	plist->set ("testRankRevealing", testRankRevealing, 
		    "Whether to test TSQR's rank-revealing capability.");
	plist->set ("printFieldNames", printFieldNames, 
		    "Whether to print field names (this is only done once, "
		    "for all Scalar types tested).");
	plist->set ("printResults", printResults, 
		    "Whether to print test results.");
	plist->set ("failIfInaccurate", failIfInaccurate,
		    "Whether to fail the test if the factorization "
		    "is not sufficiently accurate.");
	plist->set ("debug", debug, 
		    "Whether to print debugging output.");
	return plist;
      }

      /// \brief Run TsqrVerifier<T>::run() for every type in the type list.
      ///
      /// TypeListType should be either a \c NullCons (representing an
      /// empty type list, in which case this function does nothing),
      /// or a \c Cons (whose CarType is a Scalar type to test, and
      /// whose CdrType is either a NullCons or a Cons).
      /// 
      /// \param testParams [in/out] List of parameters for all tests
      ///   to run.  Call \c getValidParameterList() to get a valid
      ///   list of parameters with default values and documentation.
      ///
      template<class TypeListType>
      void 
      run (const Teuchos::RCP<Teuchos::ParameterList>& testParams)
      {
	// Using a class with a static method is a way to implement
	// "partial specialization of function templates" (which by
	// itself is not allowed in C++).
	typedef FullTsqrVerifierCallerImpl<TypeListType> impl_type;
	impl_type::run (comm_, node_, testParams, randomSeed_);
      }


      /// \brief Full constructor.
      ///
      /// \param comm [in] Communicator (with one or more processes)
      ///   over which to perform tests.
      ///
      /// \param node [in] Kokkos Node instance.
      ///
      /// \param randomSeed [in] The seed for LAPACK's pseudorandom
      ///   number generator.  An array of four integers, satisfying
      ///   the requirements of LAPACK's _LARNV routines.  The array
      ///   elements must be in [0,4095], and the last element
      ///   (iseed[3]) must be odd.  Call \c defaultRandomSeed() for a
      ///   constant default value (if you want the same results each
      ///   time; not "random" but reproducible).
      FullTsqrVerifierCaller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
			      const Teuchos::RCP<const node_type>& node,
			      const std::vector<int>& randomSeed) :
	comm_ (comm), node_ (node), randomSeed_ (validateRandomSeed (randomSeed))
      {}

      /// \brief One-argument constructor.
      ///
      /// Fills in defaults for the other arguments that the full
      /// constructor would take.
      ///
      /// \param comm [in] Communicator (with one or more processes)
      ///   over which to perform tests.
      FullTsqrVerifierCaller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm) :
	comm_ (comm),
	node_ (getNode<node_type> (getValidNodeParameters<node_type> ())),
	randomSeed_ (defaultRandomSeed ())
      {}

      /// \brief Two-argument constructor.
      ///
      /// Fills in defaults for the other arguments that the full
      /// constructor would take.
      ///
      /// \param comm [in] Communicator (with one or more processes)
      ///   over which to perform tests.
      ///
      /// \param node [in] Kokkos Node instance.
      FullTsqrVerifierCaller (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
			      const Teuchos::RCP<const node_type>& node) :
	comm_ (comm),
	node_ (node),
	randomSeed_ (defaultRandomSeed ())
      {}

      //! Validate the given random seed.
      static std::vector<int> 
      validateRandomSeed (const std::vector<int>& seed) 
      {
	TEST_FOR_EXCEPTION(seed.size() < 4, std::invalid_argument,
			   "Invalid random seed: Need an array of four integers.");
	for (std::vector<int>::size_type k = 0; k < seed.size(); ++k)
	  {
	    TEST_FOR_EXCEPTION(seed[k] < 0 || seed[k] > 4095,
			       std::invalid_argument,
			       "Invalid random seed: Each of the four integers must be in [0, 4095].");
	  }
	TEST_FOR_EXCEPTION(seed[3] % 2 != 1, std::invalid_argument,
			   "Invalid random seed: The last of the four integers must be odd.");
	return seed;
      }

      //! Default random seed.
      static std::vector<int> 
      defaultRandomSeed () 
      {
	std::vector<int> seed (4);
	seed[0] = 0;
	seed[1] = 0;
	seed[2] = 0;
	seed[3] = 1;
	return seed;
      }

    private:
      /// \brief Communicator over which to perform tests.
      ///
      /// This communicator may include one or more processes.
      /// MPI is not required (it may be a "serial communicator").
      Teuchos::RCP<const Teuchos::Comm<int> > comm_;

      //! Kokkos Node instance.
      Teuchos::RCP<const node_type> node_;

      /// \brief The seed for LAPACK's pseudorandom number generator.
      /// 
      /// Array of four integers, satisfying the requirements of
      /// LAPACK's _LARNV routines.  The array elements must be in
      /// [0,4095], and the last element (iseed[3]) must be odd.
      std::vector<int> randomSeed_;
    };

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_FullTsqrTest_hpp

