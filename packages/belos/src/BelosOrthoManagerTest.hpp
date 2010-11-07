//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2010 Sandia Corporation
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

/// \file BelosOrthoManagerTest.hpp
/// \brief Tests for Belos::OrthoManager and Belos::MatOrthoManager subclasses
///

#include "BelosConfigDefs.hpp"
#include "BelosOutputManager.hpp"
#include "BelosOrthoManager.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <stdexcept>

using std::endl;

namespace Belos {
  namespace Test {

    /// \class OrthoManagerTester
    /// \brief Wrapper around OrthoManager test functionality
    ///
    template< class Scalar, class MV >
    class OrthoManagerTester {
    private:
      typedef typename Teuchos::Array< Teuchos::RCP< MV > >::size_type size_type;

    public:
      typedef Scalar scalar_type;
      typedef Teuchos::ScalarTraits< scalar_type > SCT;
      typedef typename SCT::magnitudeType magnitude_type;
      typedef Teuchos::ScalarTraits< magnitude_type > SMT;
      typedef Belos::MultiVecTraits< scalar_type, MV > MVT;
      typedef Teuchos::SerialDenseMatrix< int, scalar_type > serial_matrix_type;

      /// \brief Run all the tests
      ///
      /// \param OM [in/out] OrthoManager subclass instance to test
      /// \param isRankRevealing [in] Whether that OrthoManager 
      ///   subclass instance has a true rank-revealing capability.
      ///   If not, we do not test it on rank-deficient vectors.
      /// \param S [in/out] Multivector instance
      /// \param sizeX1 [in] Number of columns in X1 (a multivector
      ///   instance created internally for tests)
      /// \param sizeX2 [in] Number of columns in X2 (a multivector
      ///   instance created internally for tests)
      /// \param MyOM [out] Output manager for handling local output.
      ///   In Anasazi, this class is called BasicOutputManager.  In
      ///   Belos, this class is called OutputManager.
      ///
      /// \return Number of tests that failed (zero means success)
      static int
      runTests (const Teuchos::RCP< Belos::OrthoManager< Scalar, MV > >& OM,
		const bool isRankRevealing,
		const Teuchos::RCP< MV >& S,
		const int sizeX1,
		const int sizeX2,
		const Teuchos::RCP< Belos::OutputManager< Scalar > >& MyOM)
      {
	using Teuchos::Array;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::tuple;

	// Number of tests that have failed thus far.
	int numFailed = 0; 

	// Relative tolerance against which all tests are performed.
	const magnitude_type TOL = 1.0e-12;
	// Absolute tolerance constant
	//const magnitude_type ATOL = 10;

	const scalar_type ZERO = SCT::zero();
	const scalar_type ONE = SCT::one();

	// Debug output stream
	std::ostream& debugOut = MyOM->stream(Debug);

	// Number of columns in the input "prototype" multivector S.
	const int sizeS = MVT::GetNumberVecs (*S);

	// Create multivectors X1 and X2, using the same map as multivector
	// S.  Then, test orthogonalizing X2 against X1.  After doing so, X1
	// and X2 should each be M-orthonormal, and should be mutually
	// M-orthogonal.
	debugOut << "Generating X1,X2 for testing... ";
	RCP< MV > X1 = MVT::Clone (*S, sizeX1);
	RCP< MV > X2 = MVT::Clone (*S, sizeX2);
	debugOut << "done." << endl;
	{
	  magnitude_type err;

	  //
	  // Fill X1 with random values, and test the normalization error.
	  //
	  debugOut << "Filling X2 with random values... ";
	  MVT::MvRandom(*X1);
	  debugOut << "done." << endl
		   << "Calling normalize() on X1... ";
	  // The Anasazi and Belos OrthoManager interfaces differ.
	  // For example, Anasazi's normalize() method accepts either
	  // one or two arguments, whereas Belos' normalize() requires
	  // two arguments.
	  const int initialX1Rank = OM->normalize(*X1, Teuchos::null);
	  TEST_FOR_EXCEPTION(initialX1Rank != sizeX1, 
			     std::runtime_error, 
			     "normalize(X1) returned rank "
			     << initialX1Rank << " from " << sizeX1
			     << " vectors. Cannot continue.");
	  debugOut << "done." << endl 
		   << "Calling orthonormError() on X1... ";
	  err = OM->orthonormError(*X1);
	  TEST_FOR_EXCEPTION(err > TOL, std::runtime_error,
			     "After normalize(X1), orthonormError(X1) = " 
			     << err << " > TOL = " << TOL);
	  debugOut << "done: ||<X1,X1> - I|| = " << err << endl;

	  //
	  // Fill X2 with random values, project against X1 and normalize,
	  // and test the orthogonalization error.
	  //
	  debugOut << "Filling X1 with random values... ";
	  MVT::MvRandom(*X2);
	  debugOut << "done." << endl
		   << "Calling projectAndNormalize(X2,X1)... ";
	  // The projectAndNormalize() interface also differs between 
	  // Anasazi and Belos.  Anasazi's projectAndNormalize() puts 
	  // the multivector and the array of multivectors first, and
	  // the (array of) SerialDenseMatrix arguments (which are 
	  // optional) afterwards.  Belos puts the (array of) 
	  // SerialDenseMatrix arguments in the middle, and they are 
	  // not optional.
	  int initialX2Rank;
	  {
	    // FIXME (mfh 04 Nov 2010) TsqrOrthoManager resizes C
	    // automatically.  The other OrthoManager subclasses in
	    // Belos don't seem to do this.  The commented-out line
	    // works fine with TsqrOrthoManager, but not with the
	    // others.  The current line also doesn't work, but for a
	    // different reason: the OrthoManager tries to access C[0]
	    // and fails due to a Teuchos::NullReferenceError.
	    //Array< RCP< serial_matrix_type > > C;
	    Array< RCP< serial_matrix_type > > C (1);
	    RCP< serial_matrix_type > B = Teuchos::null;
	    initialX2Rank = 
	      OM->projectAndNormalize (*X2, C, B, tuple< RCP< const MV > >(X1));
	  }
	  TEST_FOR_EXCEPTION(initialX2Rank != sizeX2, 
			     std::runtime_error, 
			     "projectAndNormalize(X2,X1) returned rank " 
			     << initialX2Rank << " from " << sizeX2 
			     << " vectors. Cannot continue.");
	  debugOut << "done." << endl
		   << "Calling orthonormError() on X2... ";
	  err = OM->orthonormError (*X2);
	  TEST_FOR_EXCEPTION(err > TOL,
			     std::runtime_error,
			     "projectAndNormalize(X2,X1) did not meet tolerance: "
			     "orthonormError(X2) = " << err << " > TOL = " << TOL);
	  debugOut << "done: || <X2,X2> - I || = " << err << endl
		   << "Calling orthogError(X2, X1)... ";
	  err = OM->orthogError (*X2, *X1);
	  TEST_FOR_EXCEPTION(err > TOL,
			     std::runtime_error,
			     "projectAndNormalize(X2,X1) did not meet tolerance: "
			     "orthogError(X2,X1) = " << err << " > TOL = " << TOL);
	  debugOut << "done: || <X2,X1> || = " << err << endl;
	}

	{
	  //
	  // Test project() on a random multivector S, by projecting S
	  // against various combinations of X1 and X2.
	  //
	  MVT::MvRandom(*S);

	  debugOut << "Testing project() by projecting a random multivector S "
	    "against various combinations of X1 and X2 " << endl;
	  const int thisNumFailed = testProject(OM,S,X1,X2,MyOM);
	  numFailed += thisNumFailed;
	  if (thisNumFailed > 0)
	    debugOut << "  *** " << thisNumFailed 
		     << (thisNumFailed > 1 ? " tests" : " test") 
		     << " failed." << endl;
	}

	if (isRankRevealing)
	  {
	    // run a X1,Y2 range multivector against P_{X1,X1} P_{Y2,Y2}
	    // note, this is allowed under the restrictions on project(), 
	    // because <X1,Y2> = 0
	    // also, <Y2,Y2> = I, but <X1,X1> != I, so biOrtho must be set to false
	    // it should require randomization, as 
	    // P_{X1,X1} P_{Y2,Y2} (X1*C1 + Y2*C2) = P_{X1,X1} X1*C1 = 0
	    serial_matrix_type C1(sizeX1,sizeS), C2(sizeX2,sizeS);
	    C1.random();
	    C2.random();
	    MVT::MvTimesMatAddMv(ONE,*X1,C1,ZERO,*S);
	    MVT::MvTimesMatAddMv(ONE,*X2,C2,ONE,*S);

	    debugOut << "Testing project() by projecting [X1 X2]-range multivector "
	      "against P_X1 P_X2 " << endl;
	    const int thisNumFailed = testProject(OM,S,X1,X2,MyOM);
	    numFailed += thisNumFailed;
	    if (thisNumFailed > 0)
	      debugOut << "  *** " << thisNumFailed 
		       << (thisNumFailed > 1 ? " tests" : " test") 
		       << " failed." << endl;
	  }

	// This test is only distinct from the rank-1 multivector test
	// (below) if S has at least 3 columns.
	if (isRankRevealing && sizeS > 2) 
	  {
	    MVT::MvRandom(*S);
	    RCP<MV> mid = MVT::Clone(*S,1);
	    serial_matrix_type c(sizeS,1);
	    MVT::MvTimesMatAddMv(ONE,*S,c,ZERO,*mid);
	    std::vector<int> ind(1); 
	    ind[0] = sizeS-1;
	    MVT::SetBlock(*mid,ind,*S);

	    debugOut << "Testing normalize() on a rank-deficient multivector " << endl;
	    const int thisNumFailed = testNormalize(OM,S,MyOM);
	    numFailed += thisNumFailed;
	    if (thisNumFailed > 0)
	      debugOut << "  *** " << thisNumFailed 
		       << (thisNumFailed > 1 ? " tests" : " test") 
		       << " failed." << endl;
	  }

	// This test will only exercise rank deficiency if S has at least 2
	// columns.
	if (isRankRevealing && sizeS > 1) 
	  {
	    // rank-1
	    RCP<MV> one = MVT::Clone(*S,1);
	    MVT::MvRandom(*one);
	    // put multiple of column 0 in columns 0:sizeS-1
	    for (int i=0; i<sizeS; i++) 
	      {
		std::vector<int> ind(1); 
		ind[0] = i;
		RCP<MV> Si = MVT::CloneViewNonConst(*S,ind);
		MVT::MvAddMv(SCT::random(),*one,ZERO,*one,*Si);
	      }
	    debugOut << "Testing normalize() on a rank-1 multivector " << endl;
	    const int thisNumFailed = testNormalize(OM,S,MyOM);
	    numFailed += thisNumFailed;
	    if (thisNumFailed > 0)
	      debugOut << "  *** " << thisNumFailed 
		       << (thisNumFailed > 1 ? " tests" : " test") 
		       << " failed." << endl;
	  }

	{
	  std::vector<int> ind(1); 
	  MVT::MvRandom(*S);

	  debugOut << "Testing projectAndNormalize() on a random multivector " << endl;
	  const int thisNumFailed = testProjectAndNormalize(OM,S,X1,X2,MyOM);
	  numFailed += thisNumFailed;
	  if (thisNumFailed > 0)
	    debugOut << "  *** " << thisNumFailed 
		     << (thisNumFailed > 1 ? " tests" : " test") 
		     << " failed." << endl;
	}

	if (isRankRevealing)
	  {
	    // run a X1,X2 range multivector against P_X1 P_X2
	    // this is allowed as <X1,X2> == 0
	    // it should require randomization, as 
	    // P_X1 P_X2 (X1*C1 + X2*C2) = P_X1 X1*C1 = 0
	    // and 
	    // P_X2 P_X1 (X2*C2 + X1*C1) = P_X2 X2*C2 = 0
	    serial_matrix_type C1(sizeX1,sizeS), C2(sizeX2,sizeS);
	    C1.random();
	    C2.random();
	    MVT::MvTimesMatAddMv(ONE,*X1,C1,ZERO,*S);
	    MVT::MvTimesMatAddMv(ONE,*X2,C2,ONE,*S);

	    debugOut << "Testing projectAndNormalize() by projecting [X1 X2]-range "
	      "multivector against P_X1 P_X2 " << endl;
	    const int thisNumFailed = testProjectAndNormalize(OM,S,X1,X2,MyOM);
	    numFailed += thisNumFailed;
	    if (thisNumFailed > 0)
	      debugOut << "  *** " << thisNumFailed 
		       << (thisNumFailed > 1 ? " tests" : " test") 
		       << " failed." << endl;
	  }

	// This test is only distinct from the rank-1 multivector test
	// (below) if S has at least 3 columns.
	if (isRankRevealing && sizeS > 2) 
	  {
	    MVT::MvRandom(*S);
	    RCP<MV> mid = MVT::Clone(*S,1);
	    serial_matrix_type c(sizeS,1);
	    MVT::MvTimesMatAddMv(ONE,*S,c,ZERO,*mid);
	    std::vector<int> ind(1); 
	    ind[0] = sizeS-1;
	    MVT::SetBlock(*mid,ind,*S);

	    debugOut << "Testing projectAndNormalize() on a rank-deficient "
	      "multivector " << endl;
	    const int thisNumFailed = testProjectAndNormalize(OM,S,X1,X2,MyOM);
	    numFailed += thisNumFailed;
	    if (thisNumFailed > 0)
	      debugOut << "  *** " << thisNumFailed 
		       << (thisNumFailed > 1 ? " tests" : " test") 
		       << " failed." << endl;
	  }

	// This test will only exercise rank deficiency if S has at least 2
	// columns.
	if (isRankRevealing && sizeS > 1) 
	  {
	    // rank-1
	    RCP<MV> one = MVT::Clone(*S,1);
	    MVT::MvRandom(*one);
	    // Put a multiple of column 0 in columns 0:sizeS-1.
	    for (int i=0; i<sizeS; i++) 
	      {
		std::vector<int> ind(1); 
		ind[0] = i;
		RCP<MV> Si = MVT::CloneViewNonConst(*S,ind);
		MVT::MvAddMv(SCT::random(),*one,ZERO,*one,*Si);
	      }
	    debugOut << "Testing projectAndNormalize() on a rank-1 multivector " << endl;
	    const int thisNumFailed = testProjectAndNormalize(OM,S,X1,X2,MyOM);
	    numFailed += thisNumFailed;
	    if (thisNumFailed > 0)
	      debugOut << "  *** " << thisNumFailed 
		       << (thisNumFailed > 1 ? " tests" : " test") 
		       << " failed." << endl;
	  }

	if (numFailed != 0)
	  MyOM->stream(Errors) << numFailed << " total test failures." << endl;

	return numFailed;
      }

    private:

      /// \fn MVDiff
      ///
      /// Compute and return $\|X - Y\|_F$, the Frobenius (sum of
      /// squares) norm of the difference between X and Y.
      static magnitude_type
      MVDiff (const MV& X, const MV& Y)
      {
	using Teuchos::RCP;

	const scalar_type ONE = SCT::one();
	const int numCols = MVT::GetNumberVecs(X);
	TEST_FOR_EXCEPTION( (MVT::GetNumberVecs(Y) != numCols),
			    std::logic_error,
			    "MVDiff: X and Y should have the same number of columns."
			    "  X has " << numCols << " column(s) and Y has " 
			    << MVT::GetNumberVecs(Y) << " columns." );
	// Resid := X
	RCP< MV > Resid = MVT::CloneCopy(X);
	// Resid := Resid - Y
	MVT::MvAddMv (-ONE, Y, ONE, *Resid, *Resid);

	return frobeniusNorm (*Resid);
      }


      /// \fn frobeniusNorm
      ///
      /// Compute and return the Frobenius norm of X.
      static magnitude_type
      frobeniusNorm (const MV& X)
      {
	const scalar_type ONE = SCT::one();
	const int numCols = MVT::GetNumberVecs(X);
	serial_matrix_type C (numCols, numCols);

	// $C := X^* X$
	MVT::MvTransMv (ONE, X, X, C);

	magnitude_type err (0);
	for (int i = 0; i < numCols; ++i)
	  err += SCT::magnitude (C(i,i));

	return SCT::magnitude (SCT::squareroot (err));
      }


      /// Test OrthoManager::projectAndNormalize() for the specific
      /// OrthoManager instance.
      ///
      /// \return Count of errors (should be zero)
      static int 
      testProjectAndNormalize (const Teuchos::RCP< Belos::OrthoManager< Scalar, MV > >& OM,
			       const Teuchos::RCP< const MV >& S, 
			       const Teuchos::RCP< const MV >& X1, 
			       const Teuchos::RCP< const MV >& X2,
			       const Teuchos::RCP< Belos::OutputManager< Scalar > >& MyOM)
      {
	using Teuchos::Array;
	using Teuchos::null;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::tuple;

	const scalar_type ONE = SCT::one();
	const magnitude_type ZERO = SCT::magnitude(SCT::zero());

	// Relative tolerance against which all tests are performed.
	const magnitude_type TOL = 1.0e-12;
	// Absolute tolerance constant
	const magnitude_type ATOL = 10;

	const int sizeS = MVT::GetNumberVecs(*S);
	const int sizeX1 = MVT::GetNumberVecs(*X1);
	const int sizeX2 = MVT::GetNumberVecs(*X2);
	int numerr = 0;
	std::ostringstream sout;

	//
	// output tests:
	//   <S_out,S_out> = I
	//   <S_out,X1> = 0
	//   <S_out,X2> = 0
	//   S_in = S_out B + X1 C1 + X2 C2
	// 
	// we will loop over an integer specifying the test combinations
	// the bit pattern for the different tests is listed in parenthesis
	//
	// for the projectors, test the following combinations:
	// none              (00)
	// P_X1              (01)
	// P_X2              (10)
	// P_X1 P_X2         (11)
	// P_X2 P_X1         (11)
	// the latter two should be tested to give the same answer
	//
	// for each of these, we should test with C1, C2 and B
	//
	// if hasM:
	// with and without MX1   (1--) 
	// with and without MX2  (1---) 
	// with and without MS  (1----) 
	//
	// as hasM controls the upper level bits, we need only run test cases 0-3 if hasM==false
	// otherwise, we run test cases 0-31
	//

	int numtests = 4;

	// test ortho error before orthonormalizing
	if (X1 != null) {
	  magnitude_type err = OM->orthogError(*S,*X1);
	  sout << "   || <S,X1> || before     : " << err << endl;
	}
	if (X2 != null) {
	  magnitude_type err = OM->orthogError(*S,*X2);
	  sout << "   || <S,X2> || before     : " << err << endl;
	}

	for (int t=0; t<numtests; t++) {

	  Array< RCP< const MV > > theX;
	  RCP<serial_matrix_type > B = rcp( new serial_matrix_type(sizeS,sizeS) );
	  Array<RCP<serial_matrix_type > > C;
	  if ( (t && 3) == 0 ) {
	    // neither <X1,Y1> nor <X2,Y2>
	    // C, theX and theY are already empty
	  }
	  else if ( (t && 3) == 1 ) {
	    // X1
	    theX = tuple(X1);
	    C = tuple( rcp(new serial_matrix_type(sizeX1,sizeS)) );
	  }
	  else if ( (t && 3) == 2 ) {
	    // X2
	    theX = tuple(X2);
	    C = tuple( rcp(new serial_matrix_type(sizeX2,sizeS)) );
	  }
	  else {
	    // X1 and X2, and the reverse.
	    theX = tuple(X1,X2);
	    C = tuple( rcp(new serial_matrix_type(sizeX1,sizeS)), 
		       rcp(new serial_matrix_type(sizeX2,sizeS)) );
	  }

	  // We wrap up all the OrthoManager calls in a try-catch
	  // block, in order to check whether any of the methods throw
	  // an exception.  For the tests we perform, every thrown
	  // exception is a failure.
	  try {
	    // call routine
	    // if (t && 3) == 3, {
	    //    call with reversed input: X2 X1
	    // }
	    // test all outputs for correctness
	    // test all outputs for equivalence

	    // here is where the outputs go
	    Array<RCP<MV> > S_outs;
	    Array<Array<RCP<serial_matrix_type > > > C_outs;
	    Array<RCP<serial_matrix_type > > B_outs;
	    RCP<MV> Scopy;
	    Array<int> ret_out;

	    // copies of S,MS
	    Scopy = MVT::CloneCopy(*S);
	    // randomize this data, it should be overwritten
	    B->random();
	    for (size_type i=0; i<C.size(); i++) {
	      C[i]->random();
	    }
	    // Run test.  Since S was specified by the caller and
	    // Scopy is a copy of S, we don't know what rank to expect
	    // here -- though we do require that S have rank at least
	    // one.
	    //
	    // Note that Anasazi and Belos differ, among other places, 
	    // in the order of arguments to projectAndNormalize().
	    int ret = OM->projectAndNormalize(*Scopy,C,B,theX);
	    sout << "projectAndNormalize() returned rank " << ret << endl;
	    if (ret == 0) {
	      sout << "  *** Error: returned rank is zero, cannot continue tests" << endl;
	      numerr++;
	      break;
	    }
	    ret_out.push_back(ret);
	    // projectAndNormalize() is only required to return a 
	    // basis of rank "ret"
	    // this is what we will test:
	    //   the first "ret" columns in Scopy
	    //   the first "ret" rows in B
	    // save just the parts that we want
	    // we allocate S and MS for each test, so we can save these as views
	    // however, save copies of the C and B
	    if (ret < sizeS) {
	      std::vector<int> ind(ret);
	      for (int i=0; i<ret; i++) {
		ind[i] = i;
	      }
	      S_outs.push_back( MVT::CloneViewNonConst(*Scopy,ind) );
	      B_outs.push_back( rcp( new serial_matrix_type(Teuchos::Copy,*B,ret,sizeS) ) );
	    }
	    else {
	      S_outs.push_back( Scopy );
	      B_outs.push_back( rcp( new serial_matrix_type(*B) ) );
	    }
	    C_outs.push_back( Array<RCP<serial_matrix_type > >(0) );
	    if (C.size() > 0) {
	      C_outs.back().push_back( rcp( new serial_matrix_type(*C[0]) ) );
	    }
	    if (C.size() > 1) {
	      C_outs.back().push_back( rcp( new serial_matrix_type(*C[1]) ) );
	    }

	    // do we run the reversed input?
	    if ( (t && 3) == 3 ) {
	      // copies of S,MS
	      Scopy = MVT::CloneCopy(*S);

	      // Fill the B and C[i] matrices with random data.  The
	      // data will be overwritten by projectAndNormalize().
	      // Filling these matrices here is only to catch some
	      // bugs in projectAndNormalize().
	      B->random();
	      for (size_type i=0; i<C.size(); i++) {
		C[i]->random();
	      }
	      // flip the inputs
	      theX = tuple( theX[1], theX[0] );
	      // Run test.
	      // Note that Anasazi and Belos differ, among other places, 
	      // in the order of arguments to projectAndNormalize().
	      ret = OM->projectAndNormalize(*Scopy,C,B,theX);
	      sout << "projectAndNormalize() returned rank " << ret << endl;
	      if (ret == 0) {
		sout << "  *** Error: returned rank is zero, cannot continue tests" << endl;
		numerr++;
		break;
	      }
	      ret_out.push_back(ret);
	      // projectAndNormalize() is only required to return a 
	      // basis of rank "ret"
	      // this is what we will test:
	      //   the first "ret" columns in Scopy
	      //   the first "ret" rows in B
	      // save just the parts that we want
	      // we allocate S and MS for each test, so we can save these as views
	      // however, save copies of the C and B
	      if (ret < sizeS) {
		std::vector<int> ind(ret);
		for (int i=0; i<ret; i++) {
		  ind[i] = i;
		}
		S_outs.push_back( MVT::CloneViewNonConst(*Scopy,ind) );
		B_outs.push_back( rcp( new serial_matrix_type(Teuchos::Copy,*B,ret,sizeS) ) );
	      }
	      else {
		S_outs.push_back( Scopy );
		B_outs.push_back( rcp( new serial_matrix_type(*B) ) );
	      }
	      C_outs.push_back( Array<RCP<serial_matrix_type > >() );
	      // reverse the Cs to compensate for the reverse projectors
	      C_outs.back().push_back( rcp( new serial_matrix_type(*C[1]) ) );
	      C_outs.back().push_back( rcp( new serial_matrix_type(*C[0]) ) );
	      // flip the inputs back
	      theX = tuple( theX[1], theX[0] );
	    }


	    // test all outputs for correctness
	    for (size_type o=0; o<S_outs.size(); o++) {
	      // S^T M S == I
	      {
		magnitude_type err = OM->orthonormError(*S_outs[o]);
		if (err > TOL) {
		  sout << endl
		       << "  *** Test (number " << (t+1) << " of " << numtests 
		       << " total tests) failed: Tolerance exceeded!  Error = "
		       << err << " > TOL = " << TOL << "." 
		       << endl << endl;
		  numerr++;
		}
		sout << "   || <S,S> - I || after  : " << err << endl;
	      }
	      // S_in = X1*C1 + C2*C2 + S_out*B
	      {
		RCP<MV> tmp = MVT::Clone(*S,sizeS);
		MVT::MvTimesMatAddMv(ONE,*S_outs[o],*B_outs[o],ZERO,*tmp);
		if (C_outs[o].size() > 0) {
		  MVT::MvTimesMatAddMv(ONE,*X1,*C_outs[o][0],ONE,*tmp);
		  if (C_outs[o].size() > 1) {
		    MVT::MvTimesMatAddMv(ONE,*X2,*C_outs[o][1],ONE,*tmp);
		  }
		}
		magnitude_type err = MVDiff(*tmp,*S);
		if (err > ATOL*TOL) {
		  sout << endl
		       << "  *** Test (number " << (t+1) << " of " << numtests 
		       << " total tests) failed: Tolerance exceeded!  Error = "
		       << err << " > ATOL*TOL = " << (ATOL*TOL) << "."
		       << endl << endl;
		  numerr++;
		}
		sout << "  " << t << "|| S_in - X1*C1 - X2*C2 - S_out*B || : " << err << endl;
	      }
	      // <X1,S> == 0
	      if (theX.size() > 0 && theX[0] != null) {
		magnitude_type err = OM->orthogError(*theX[0],*S_outs[o]);
		if (err > TOL) {
		  sout << endl
		       << "  *** Test (number " << (t+1) << " of " << numtests 
		       << " total tests) failed: Tolerance exceeded!  Error = "
		       << err << " > TOL = " << TOL << "."
		       << endl << endl;
		  numerr++;
		}
		sout << "  " << t << "|| <X[0],S> || after      : " << err << endl;
	      }
	      // <X2,S> == 0
	      if (theX.size() > 1 && theX[1] != null) {
		magnitude_type err = OM->orthogError(*theX[1],*S_outs[o]);
		if (err > TOL) {
		  sout << endl
		       << "  *** Test (number " << (t+1) << " of " << numtests 
		       << " total tests) failed: Tolerance exceeded!  Error = "
		       << err << " > TOL = " << TOL << "."
		       << endl << endl;
		  numerr++;
		}
		sout << "  " << t << "|| <X[1],S> || after      : " << err << endl;
	      }
	    }
	  }
	  catch (Belos::OrthoError& e) {
	    sout << "  *** Error: OrthoManager threw exception: " << e.what() << endl;
	    numerr++;
	  }

	} // test for

	// FIXME (mfh 05 Nov 2010) Since Belos::MsgType is an enum,
	// doing bitwise logical computations on Belos::MsgType values
	// (such as "Debug | Errors") and passing the result into
	// MyOM->stream() confuses the compiler.  As a result, we have
	// to do weird casts to make it work.
	const int msgType = (numerr > 0) ? 
	  (static_cast<int>(Debug) | static_cast<int>(Errors)) :
	  static_cast<int>(Debug);

	// We report debug-level messages always.  We also report
	// errors if at least one test failed.
	MyOM->stream(static_cast< MsgType >(msgType)) << sout.str() << endl;
	return numerr;
      }

      /// Test OrthoManager::normalize() for the specific OrthoManager
      /// instance.
      ///
      /// \return Count of errors (should be zero)
      static int 
      testNormalize (const Teuchos::RCP< Belos::OrthoManager< Scalar, MV > >& OM, 
		     const Teuchos::RCP< const MV >& S,
		     const Teuchos::RCP< Belos::OutputManager< Scalar > >& MyOM)
      {
	using Teuchos::Array;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::tuple;

	const scalar_type ONE = SCT::one();
	std::ostringstream sout;
	// Total number of failed tests in this call of this routine.
	int numerr = 0;

	// Relative tolerance against which all tests are performed.
	// We are measuring things in the Frobenius norm $\| \cdot \|_F$.
	// The following bounds hold for all $m \times n$ matrices $A$:
	// \[
	// \|A\|_2 \leq \|A\|_F \leq \sqrt{r} \|A\|_2,
	// \]
	// where $r$ is the (column) rank of $A$.  We bound this above
	// by the number of columns in $A$.  
	//
	// An accurate normalization in the Euclidean norm of a matrix
	// $A$ with at least as many rows m as columns n, should
	// produce orthogonality $\|Q^* Q - I\|_2$ less than a factor
	// of machine precision times a low-order polynomial in m and
	// n, and residual $\|A - Q B\|_2$ (where $A = Q B$ is the
	// computed normalization) less than that bound times the norm
	// of $A$.
	//
	// Since we are measuring both of these quantitites in the
	// Frobenius norm instead, we should scale this bound by
	// $\sqrt{n}$.

	const int numRows = MVT::GetVecLength(*S);
	const int numCols = MVT::GetNumberVecs(*S);
	const int sizeS = MVT::GetNumberVecs(*S);

	// A good heuristic is to scale the bound by the square root
	// of the number of floating-point operations.  One could
	// perhaps support this theoretically, since we are using
	// uniform random test problems.
	const magnitude_type fudgeFactor = 
	  SMT::squareroot(magnitude_type(numRows) * 
			  magnitude_type(numCols) * 
			  magnitude_type(numCols));
	const magnitude_type TOL = SMT::eps() * fudgeFactor *
	  SMT::squareroot(magnitude_type(numCols));
	// Absolute tolerance scaling: the Frobenius norm of the test
	// matrix S.  TOL*ATOL is the absolute tolerance for the
	// residual $\|A - Q*B\|_F$.
	
	std::cerr << "+++ COMPUTING FROB NORM +++" << endl;
	const magnitude_type ATOL = frobeniusNorm (*S);
	std::cerr << "+++ DONE COMPUTING FROB NORM +++" << endl;

	sout << "The test matrix S has Frobenius norm " << ATOL 
	     << ", and the relative error tolerance is TOL = " 
	     << TOL << "." << endl;

	const int numtests = 1;
	for (int t = 0; t < numtests; ++t) {

	  try {
	    // call routine
	    // test all outputs for correctness

	    // S_copy gets a copy of S; we normalize in place, so we
	    // need a copy to check whether the normalization
	    // succeeded.
	    RCP< MV > S_copy = MVT::CloneCopy (*S);

	    // Matrix of coefficients from the normalization.
	    RCP< serial_matrix_type > B (new serial_matrix_type (sizeS, sizeS));
	    // The contents of B will be overwritten, but fill with
	    // random data just to make sure that the normalization
	    // operated on all the elements of B on which it should
	    // operate.
	    B->random();

	    const int reportedRank = OM->normalize (*S_copy, B);
	    sout << "normalize() returned rank " << reportedRank << endl;
	    if (reportedRank == 0) {
	      sout << "  *** Error: Cannot continue, since normalize() "
		"reports that S has rank 0" << endl;
	      numerr++;
	      break;
	    }
	    //
	    // We don't know in this routine whether the input
	    // multivector S has full rank; it is only required to
	    // have nonzero rank.  Thus, we extract the first
	    // reportedRank columns of S_copy and the first
	    // reportedRank rows of B, and perform tests on them.
	    //

	    // Construct S_view, a view of the first reportedRank
	    // columns of S_copy.
	    std::vector<int> indices (reportedRank);
	    for (int j = 0; j < reportedRank; ++j)
	      indices[j] = j;
	    RCP< MV > S_view = MVT::CloneViewNonConst (*S_copy, indices);
	    // Construct B_top, a copy of the first reportedRank rows
	    // of B.
	    //
	    // NOTE: We create this as a copy and not a view, because
	    // otherwise it would not be safe with respect to RCPs.
	    // This is because serial_matrix_type uses raw pointers
	    // inside, so that a view would become invalid when B
	    // would fall out of scope.
	    RCP< serial_matrix_type > B_top (new serial_matrix_type (Teuchos::Copy, *B, reportedRank, sizeS));

	    // Check ||<S_view,S_view> - I||
	    {
	      const magnitude_type err = OM->orthonormError(*S_view);
	      if (err > TOL) {
		sout << "  *** Error: Tolerance exceeded: err = " 
		     << err << " > TOL = " << TOL << endl;
		numerr++;
	      }
	      sout << "   || <S,S> - I || after  : " << err << endl;
	    }
	    // Check the residual ||Residual|| = ||S_view * B_top -
	    // S_orig||, where S_orig is a view of the first
	    // reportedRank columns of S.
	    {
	      // Residual is allocated with reportedRank columns.  It
	      // will contain the result of testing the residual error
	      // of the normalization (i.e., $\|S - S_in*B\|$).  It
	      // should have the dimensions of S.  Its initial value
	      // is a copy of the first reportedRank columns of S.
	      RCP< MV > Residual = MVT::CloneCopy (*S);

	      // Residual := Residual - S_view * B_view
	      MVT::MvTimesMatAddMv (-ONE, *S_view, *B_top, ONE, *Residual);

	      // Compute ||Residual||
	      const magnitude_type err = frobeniusNorm (*Residual);
	      if (err > ATOL*TOL) {
		sout << "  *** Error: Tolerance exceeded: err = " 
		     << err << " > ATOL*TOL = " << (ATOL*TOL) << endl;
		numerr++;
	      }
	      sout << "  " << t << "|| S - Q*B || : " << err << endl;
	    }
	  }
	  catch (Belos::OrthoError& e) {
	    sout << "  *** Error: the OrthoManager's normalize() method "
	      "threw an exception: " << e.what() << endl;
	    numerr++;
	  }

	} // test for

	const MsgType type = (numerr == 0) ? Debug : static_cast<MsgType> (static_cast<int>(Errors) | static_cast<int>(Debug));
	MyOM->stream(type) << sout.str();
	MyOM->stream(type) << endl;

	return numerr;
      }

      /// Test OrthoManager::project() for the specific OrthoManager instance.
      ///
      /// \return Count of errors (should be zero)
      static int 
      testProject (const Teuchos::RCP< Belos::OrthoManager< Scalar, MV > > OM, 
		   const Teuchos::RCP< const MV >& S, 
		   const Teuchos::RCP< const MV >& X1, 
		   const Teuchos::RCP< const MV >& X2,
		   const Teuchos::RCP< Belos::OutputManager< Scalar > >& MyOM)
      {
	using Teuchos::Array;
	using Teuchos::null;
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::tuple;

	const scalar_type ONE = SCT::one();

	// Relative tolerance against which all tests are performed.
	const magnitude_type TOL = 1.0e-12;
	// Absolute tolerance constant
	const magnitude_type ATOL = 10;

	const int sizeS = MVT::GetNumberVecs(*S);
	const int sizeX1 = MVT::GetNumberVecs(*X1);
	const int sizeX2 = MVT::GetNumberVecs(*X2);
	int numerr = 0;
	std::ostringstream sout;

	//
	// Output tests:
	//   <S_out,X1> = 0
	//   <S_out,X2> = 0
	//   S_in = S_out + X1 C1 + X2 C2
	// 
	// We will loop over an integer specifying the test combinations.
	// The bit pattern for the different tests is listed in parentheses.
	//
	// For the projectors, test the following combinations:
	// none              (00)
	// P_X1              (01)
	// P_X2              (10)
	// P_X1 P_X2         (11)
	// P_X2 P_X1         (11)
	// The latter two should be tested to give the same result.
	//
	// For each of these, we should test with C1 and C2:
	//
	// if hasM:
	// with and without MX1   (1--) 
	// with and without MX2  (1---) 
	// with and without MS  (1----) 
	//
	// As hasM controls the upper level bits, we need only run test
	// cases 0-3 if hasM==false.  Otherwise, we run test cases 0-31.
	//

	int numtests = 8;

	// test ortho error before orthonormalizing
	if (X1 != null) {
	  magnitude_type err = OM->orthogError(*S,*X1);
	  sout << "   || <S,X1> || before     : " << err << endl;
	}
	if (X2 != null) {
	  magnitude_type err = OM->orthogError(*S,*X2);
	  sout << "   || <S,X2> || before     : " << err << endl;
	}

	for (int t = 0; t < numtests; ++t) 
	  {
	    Array< RCP< const MV > > theX;
	    Array< RCP< serial_matrix_type > > C;
	    if ( (t && 3) == 0 ) {
	      // neither X1 nor X2
	      // C and theX are already empty
	    }
	    else if ( (t && 3) == 1 ) {
	      // X1
	      theX = tuple(X1);
	      C = tuple( rcp(new serial_matrix_type(sizeX1,sizeS)) );
	    }
	    else if ( (t && 3) == 2 ) {
	      // X2
	      theX = tuple(X2);
	      C = tuple( rcp(new serial_matrix_type(sizeX2,sizeS)) );
	    }
	    else {
	      // X1 and X2, and the reverse.
	      theX = tuple(X1,X2);
	      C = tuple( rcp(new serial_matrix_type(sizeX1,sizeS)), 
			 rcp(new serial_matrix_type(sizeX2,sizeS)) );
	    }

	    try {
	      // call routine
	      // if (t && 3) == 3, {
	      //    call with reversed input: X2 X1
	      // }
	      // test all outputs for correctness
	      // test all outputs for equivalence

	      // here is where the outputs go
	      Array< RCP< MV > > S_outs;
	      Array< Array< RCP< serial_matrix_type > > > C_outs;
	      RCP< MV > Scopy;

	      // copies of S,MS
	      Scopy = MVT::CloneCopy(*S);
	      // randomize this data, it should be overwritten
	      for (size_type i = 0; i < C.size(); ++i) {
		C[i]->random();
	      }
	      // Run test.
	      // Note that Anasazi and Belos differ, among other places, 
	      // in the order of arguments to project().
	      OM->project(*Scopy,C,theX);
	      // we allocate S and MS for each test, so we can save these as views
	      // however, save copies of the C
	      S_outs.push_back( Scopy );
	      C_outs.push_back( Array< RCP< serial_matrix_type > >(0) );
	      if (C.size() > 0) {
		C_outs.back().push_back( rcp( new serial_matrix_type(*C[0]) ) );
	      }
	      if (C.size() > 1) {
		C_outs.back().push_back( rcp( new serial_matrix_type(*C[1]) ) );
	      }

	      // do we run the reversed input?
	      if ( (t && 3) == 3 ) {
		// copies of S,MS
		Scopy = MVT::CloneCopy(*S);
		// randomize this data, it should be overwritten
		for (size_type i = 0; i < C.size(); ++i) {
		  C[i]->random();
		}
		// flip the inputs
		theX = tuple( theX[1], theX[0] );
		// Run test.
		// Note that Anasazi and Belos differ, among other places, 
		// in the order of arguments to project().
		OM->project(*Scopy,C,theX);
		// we allocate S and MS for each test, so we can save these as views
		// however, save copies of the C
		S_outs.push_back( Scopy );
		// we are in a special case: P_X1 and P_X2, so we know we applied 
		// two projectors, and therefore have two C[i]
		C_outs.push_back( Array<RCP<serial_matrix_type > >() );
		// reverse the Cs to compensate for the reverse projectors
		C_outs.back().push_back( rcp( new serial_matrix_type(*C[1]) ) );
		C_outs.back().push_back( rcp( new serial_matrix_type(*C[0]) ) );
		// flip the inputs back
		theX = tuple( theX[1], theX[0] );
	      }

	      // test all outputs for correctness
	      for (size_type o = 0; o < S_outs.size(); ++o) {
		// S_in = X1*C1 + C2*C2 + S_out
		{
		  RCP<MV> tmp = MVT::CloneCopy(*S_outs[o]);
		  if (C_outs[o].size() > 0) {
		    MVT::MvTimesMatAddMv(ONE,*X1,*C_outs[o][0],ONE,*tmp);
		    if (C_outs[o].size() > 1) {
		      MVT::MvTimesMatAddMv(ONE,*X2,*C_outs[o][1],ONE,*tmp);
		    }
		  }
		  magnitude_type err = MVDiff(*tmp,*S);
		  if (err > ATOL*TOL) {
		    sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
		    numerr++;
		  }
		  sout << "  " << t << "|| S_in - X1*C1 - X2*C2 - S_out || : " << err << endl;
		}
		// <X1,S> == 0
		if (theX.size() > 0 && theX[0] != null) {
		  magnitude_type err = OM->orthogError(*theX[0],*S_outs[o]);
		  if (err > TOL) {
		    sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
		    numerr++;
		  }
		  sout << "  " << t << "|| <X[0],S> || after      : " << err << endl;
		}
		// <X2,S> == 0
		if (theX.size() > 1 && theX[1] != null) {
		  magnitude_type err = OM->orthogError(*theX[1],*S_outs[o]);
		  if (err > TOL) {
		    sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
		    numerr++;
		  }
		  sout << "  " << t << "|| <X[1],S> || after      : " << err << endl;
		}
	      }

	      // test all outputs for equivalence
	      // check all combinations:
	      //    output 0 == output 1
	      //    output 0 == output 2
	      //    output 1 == output 2
	      for (size_type o1=0; o1<S_outs.size(); o1++) {
		for (size_type o2=o1+1; o2<S_outs.size(); o2++) {
		  // don't need to check MS_outs because we check 
		  //   S_outs and MS_outs = M*S_outs
		  // don't need to check C_outs either
		  //   
		  // check that S_outs[o1] == S_outs[o2]
		  magnitude_type err = MVDiff(*S_outs[o1],*S_outs[o2]);
		  if (err > TOL) {
		    sout << "    vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
		    numerr++;
		  }
		}
	      }

	    }
	    catch (Belos::OrthoError& e) {
	      sout << "   -------------------------------------------         project() threw exception" << endl;
	      sout << "   Error: " << e.what() << endl;
	      numerr++;
	    }
	  } // test for

	MsgType type = Debug;
	if (numerr>0) type = Errors;
	MyOM->stream(type) << sout.str();
	MyOM->stream(type) << endl;

	return numerr;
      }


    };



  } // namespace Test
} // namespace Belos


