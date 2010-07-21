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

/// \file cxx_main_tsqr.cpp 
/// \brief Test TsqrOrthoManager and TsqrMatOrthoManager
///
/// Test the OrthoManager interface to TsqrOrthoManager and
/// TsqrMatOrthoManager, using Tpetra::MultiVector as the multivector
/// implementation.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziTsqrOrthoManager.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <iohb.h>

#include <complex>
#include <stdexcept>

using namespace Anasazi;
using namespace Teuchos;

using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::Map;
using Tpetra::MultiVector;
using std::cout;
using std::endl;
using std::vector;

typedef double                              ST;
typedef ScalarTraits< ST >                 SCT;
typedef SCT::magnitudeType                  MT;
typedef MultiVector< ST, int >              MV;
typedef Operator< ST, int >                 OP;
typedef MultiVecTraits< ST, MV >           MVT;
typedef OperatorTraits< ST, MV, OP >       OPT;

typedef SerialDenseMatrix< int, ST > serial_matrix_type;

// this is the tolerance that all tests are performed against
const MT TOL = 1.0e-12;
const MT ATOL = 10;

// declare an output manager for handling local output
RCP< Anasazi::BasicOutputManager<ST> > MyOM;

////////////////////////////////////////////////////////////////////////////////
// Some forward declarations
////////////////////////////////////////////////////////////////////////////////

static int 
testProject (RCP<OrthoManager<ST,MV> > OM, 
	     RCP<const MV> S, 
	     RCP<const MV> X1, 
	     RCP<const MV> X2);

static int 
testNormalize (RCP<OrthoManager<ST,MV> > OM, 
	       RCP<const MV> S);

static int 
testProjectAndNormalize (RCP<OrthoManager<ST,MV> > OM, 
			 RCP<const MV> S, 
			 RCP<const MV> X1, 
			 RCP<const MV> X2);

/// Compute and return $\sum_{j=1}^n \| X(:,j) - Y(:,j) \|_2$, where
/// $n$ is the number of columns in X.
static MT 
MVDiff (const MV& X, 
	const MV& Y);


/// Valid command-line parameter values for the OrthoManager subclass
/// to test.
static const char* validOrthoManagers[] = {
  "Tsqr",
  "SVQB",
  "Basic",
  "ICGS"
};
/// Number of valid command-line parameter values for the OrthoManager
/// subclass to test.  Must be at least one.
static const int numValidOrthoManagers = 4;

/// Return a list (as a string) of valid command-line parameter values
/// for the OrthoManager subclass to test.
static std::string
printValidOrthoManagerList ()
{
  TEST_FOR_EXCEPTION( numValidOrthoManagers <= 0,
		      std::logic_error,
		      "Invalid number " 
		      << numValidOrthoManagers 
		      << " of OrthoManager command-line options" );
  if (numValidOrthoManagers > 1)
    {
      for (int k = 0; k < numValidOrthoManagers - 1; ++k)
	os << "\"" << validOrthoManagers[k] << "\", ";
      os << "or ";
    }
  os << "\"" << validOrthoManagers[numValidOrthoManagers-1] << "\"";
}
static std::string
defaultOrthoManagerName () { return std::string("Tsqr"); }


///
/// Instantiate and return an RCP to the specified OrthoManager
/// subclass.
///
static RCP< OrthoManager<ST,MV> >
getOrthoManager (const std::string& ortho)
{
    if (ortho == "SVQB") {
      return rcp (new SVQBOrthoManager<ST,MV,OP>(M));
    }
    else if (ortho == "Basic" || ortho == "basic") {
      return rcp (new BasicOrthoManager<ST,MV,OP>(M));
    }
    else if (ortho == "ICGS") {
      return rcp (new ICGSOrthoManager<ST,MV,OP>(M));
    }
    else if (ortho == "Tsqr" || ortho == "TSQR") {
      return rcp (new TsqrMatOrthoManager< ST, MV, OP > (M));
    }
    else {
      TEST_FOR_EXCEPTION( true, std::invalid_argument, 
			  "Invalid value for command-line parameter \"ortho\":"
			  " valid values are " << printValidOrthoManagerList() 
			  << "." );
    }
}


int 
main (int argc, char *argv[]) 
{
  const ST ONE = SCT::one();
  const MT ZERO = SCT::magnitude(SCT::zero());
  GlobalMPISession mpisess(&argc,&argv,&std::cout);

  int info = 0;
  int MyPID = 0;

  RCP< const Teuchos::Comm<int> > comm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  MyPID = rank(*comm);

  int numFailed = 0;
  bool verbose = false;
  bool debug = false;
  std::string filename;
  std::string ortho (defaultOrthoManagerName());
  int numRows = 100;
  int sizeS  = 5;
  int sizeX1 = 11; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  int sizeX2 = 13; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  bool success = true;
  try {

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption ("verbose", "quiet", &verbose,
		    "Print messages and results.");
    cmdp.setOption ("debug", "nodebug", &debug,
		    "Print debugging information.");
    cmdp.setOption ("filename", &filename,
		    "Filename for Harwell-Boeing sparse matrix, used as the M "
		    "operator for the inner product.  If not provided, no "
		    "matrix is read and the Euclidean inner product is used.");
    {
      std::ostringstream os;
      os << "OrthoManager subclass to test.  There ";
      if (numValidOrthoManagers > 1)
	os << "are " << numValidOrthoManagers << "options: ";
      else {
	os << "is " << numValidOrthoManagers << "option: ";
      }
      os << printValidOrthoManagerList() << ".";
      cmdp.setOption ("ortho", &ortho, os.str());
    }
    cmdp.setOption ("numRows", &numRows, 
		    "Controls the number of rows of the test "
		    "multivectors.  If an input matrix is given, this "
		    "parameter\'s value is ignored.");
    cmdp.setOption ("sizeS", &sizeS, "Controls the number of columns of the "
		    "input multivector.");
    cmdp.setOption ("sizeX1", &sizeX1, "Controls the number of columns of the "
		    "first basis.");
    cmdp.setOption ("sizeX2", &sizeX2, "Controls the number of columns of the "
		    "second basis.  We require for simplicity of testing (the "
		    "routines do not require it) that sizeX1 >= sizeX2.");
    if (cmdp.parse (argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) 
      return -1;
    if (debug) 
      verbose = true;
    
    MyOM = rcp( new BasicOutputManager<ST>() );
    if (verbose) {
      // output in this driver will be sent to Anasazi::Warnings
      MyOM->setVerbosity(Anasazi::Warnings);
    }
    MyOM->stream(Anasazi::Warnings) << Anasazi_Version() << endl << endl;

    RCP<Map<int> > map;
    RCP<CrsMatrix<ST,int> > M;
    if (filename != "") 
      {
	int numCols = 0;
	int nnz = -1;
	int rnnzmax = 0;
	double *dvals = NULL;
	int *colptr = NULL;
	int *rowind = NULL;

	if (MyPID == 0) 
	  {
	    // Proc 0 reads the sparse matrix (stored in Harwell-Boeing
	    // format) from the file into the tuple (numRows, numCols, nnz,
	    // colptr, rowind, dvals).  The routine allocates memory for
	    // colptr, rowind, and dvals using malloc().
	    info = readHB_newmat_double (filename.c_str(), &numRows, &numCols,
					 &nnz, &colptr, &rowind, &dvals);
	    // rnnzmax := maximum number of nonzeros per row, over all
	    // rows of the sparse matrix.
	    vector<int> rnnz (numRows, 0);
	    for (int *ri = rowind; ri < rowind + nnz; ++ri) {
	      ++rnnz[*ri-1];
	    }
	    rnnzmax = *std::max_element (rnnz.begin(),rnnz.end());
	  }

	// Proc 0 now broadcasts the sparse matrix data to the other
	// process(es).  First things broadcast are info and nnz, which
	// tell the other process(es) whether reading the sparse matrix
	// succeeded.  (info should be nonzero if so.  The
	// Harwell-Boeing routines return "C boolean true" rather than
	// the POSIX-standard "zero for success.")
	Teuchos::broadcast (*comm, 0, &info);
	Teuchos::broadcast (*comm, 0, &nnz);
	Teuchos::broadcast (*comm, 0, &numRows);
	Teuchos::broadcast (*comm, 0, &numCols);
	Teuchos::broadcast (*comm, 0, &rnnzmax);
	if (info == 0 || nnz < 0) 
	  {
	    if (MyPID == 0) {
	      cout << "Error reading Harwell-Boeing sparse matrix file \"" 
		   << filename << "\"" 
		   << endl
		   << "End Result: TEST FAILED" << endl;
	    }
	    return -1;
	  }
	else if (numRows != numcols)
	  {
	    if (MyPID == 0) {
	      cout << "Test matrix in Harwell-Boeing sparse matrix file '" 
		   << filename << "' " << "is not square: it is " << numRows 
		   << " by " << numCols
		   << endl
		   << "End Result: TEST FAILED" << endl;
	    }
	    return -1;
	  }
	else if (nnz == 0) 
	  {
	    if (MyPID == 0) {
	      cout << "Test matrix in Harwell-Boeing sparse matrix file '" 
		   << filename << "' " << "has zero nonzero values, which "
		   << "means it does not define a valid inner product." 
		   << endl
		   << "End Result: TEST FAILED" << endl;
	    }
	    return -1;
	  }
    
	// Create Tpetra::Map to represent multivectors in the range of
	// the sparse matrix.
	map = rcp (new Map<int> (numRows, 0, comm));
	M = rcp (new CrsMatrix<ST,int> (map, rnnzmax));

	if (MyPID == 0) 
	  {
	    // Convert from Harwell-Boeing format (compressed sparse
	    // column, one-indexed) to CrsMatrix format (compressed
	    // sparse row, zero-index).  We do this by iterating over
	    // all the columns of the matrix.
	    const int curNonzeroIndex = 0;
	    for (int c = 0; c < numCols; ++c) 
	      {
		for (int colnnz = 0; colnnz < colptr[c+1] - colptr[c]; ++colnnz) 
		  {
		    // Row index: *rptr - 1 (1-based -> 0-based indexing)
		    // Column index: c
		    // Value to insert there: *dptr
		    const int curGlobalRowIndex = rowind[curNonzeroIndex] - 1;
		    const ST curValue = dvals[curNonzeroIndex];
		    M->insertGlobalValues (curGlobalRowIndex, 
					   tuple (c), 
					   tuple (curValue));
		    curNonzeroIndex++;
		  }
	      }
	  }
	if (MyPID == 0) 
	  {
	    // Free memory allocated by the Harwell-Boeing input routine.
	    free (dvals);
	    dvals = NULL;
	    free (colptr);
	    colptr = NULL;
	    free (rowind);
	    rowind = NULL;
	  }
	// We're done reading in M.
	M->fillComplete();
      } // else M == null
    else {
      // Let M remain null, and allocate map using the number of rows
      // (numRows) specified on the command line.
      map = rcp (new Map<int> (numRows, 0, comm));
    }
    
    // Instantiate the specified OrthoManager subclass for testing.
    RCP< OrthoManager<ST,MV> > OM = getOrthoManager (ortho);

    // "Prototype" multivector, from which to clone other multivectors.
    RCP< MV > S = rcp (new MultiVector< ST, int > (map, sizeS));

    // Create multivectors X1 and X2, using the same map as
    // multivector S.  X1 and X2 must be M-orthonormal and mutually
    // M-orthogonal.
    MyOM->stream(Errors) << " Generating X1,X2 for testing... " << endl;
    RCP< MV > X1 = MVT::Clone (*S, sizeX1);
    RCP< MV > X2 = MVT::Clone (*S, sizeX2);
    {
      int dummy;
      MT err;

      //
      // Fill X1 with random values, and test the normalization error.
      //
      MVT::MvRandom(*X1);
      const int initialX1Rank = OM->normalize(*X1);
      TEST_FOR_EXCEPTION(initialX1Rank != sizeX1, 
			 std::runtime_error, 
			 "normalize(X1) returned rank "
			 << initialX1Rank << " from " << sizeX1
			 << " vectors. Cannot continue.");
      err = OM->orthonormError(*X1);
      TEST_FOR_EXCEPTION(err > TOL,
			 std::runtime_error,
			 "normalize(X1) did meet tolerance: "
			 "orthonormError(X1) == " << err);
      MyOM->stream(Warnings) << "   || <X1,X1> - I || : " << err << endl;

      //
      // Fill X2 with random values, project against X1 and normalize,
      // and test the orthogonalization error.
      //
      MVT::MvRandom(*X2);
      const int initialX2Rank = 
	OM->projectAndNormalize (*X2, tuple< RCP< const MV > > (X1));
      TEST_FOR_EXCEPTION(initialX2Rank != sizeX2, 
			 std::runtime_error, 
			 "projectAndNormalize(X2,X1) returned rank " 
			 << initialX2Rank << " from " << sizeX2 
			 << " vectors. Cannot continue.");
      err = OM->orthonormError(*X2);
      TEST_FOR_EXCEPTION(err > TOL,
			 std::runtime_error,
			 "projectAndNormalize(X2,X1) did not meet tolerance: "
			 "orthonormError(X2) == " << err);
      MyOM->stream(Warnings) << "   || <X2,X2> - I || : " << err << endl;
      err = OM->orthogError(*X2,*X1);
      TEST_FOR_EXCEPTION(err > TOL,
			 std::runtime_error,
			 "projectAndNormalize(X2,X1) did not meet tolerance: "
			 "orthogError(X2,X1) == " << err);
      MyOM->stream(Warnings) << "   || <X2,X1> ||     : " << err << endl;
    }
    MyOM->stream(Warnings) << endl;


    {
      //
      // Test project() on a random multivector S, by projecting S
      // against various combinations of X1 and X2.
      //
      MVT::MvRandom(*S);

      MyOM->stream(Errors) << " project(): testing on random multivector " << endl;
      numFailed += testProject(OM,S,X1,X2);
    }

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

      MyOM->stream(Errors) << " project(): testing [X1 X2]-range multivector against P_X1 P_X2 " << endl;
      numFailed += testProject(OM,S,X1,X2);
    }


    if (sizeS > 2) {
      MVT::MvRandom(*S);
      RCP<MV> mid = MVT::Clone(*S,1);
      serial_matrix_type c(sizeS,1);
      MVT::MvTimesMatAddMv(ONE,*S,c,ZERO,*mid);
      std::vector<int> ind(1); 
      ind[0] = sizeS-1;
      MVT::SetBlock(*mid,ind,*S);

      MyOM->stream(Errors) << " normalize(): testing on rank-deficient multivector " << endl;
      numFailed += testNormalize(OM,S);
    }


    if (sizeS > 1) {
      // rank-1
      RCP<MV> one = MVT::Clone(*S,1);
      MVT::MvRandom(*one);
      // put multiple of column 0 in columns 0:sizeS-1
      for (int i=0; i<sizeS; i++) {
        std::vector<int> ind(1); 
        ind[0] = i;
        RCP<MV> Si = MVT::CloneViewNonConst(*S,ind);
        MVT::MvAddMv(SCT::random(),*one,ZERO,*one,*Si);
      }

      MyOM->stream(Errors) << " normalize(): testing on rank-1 multivector " << endl;
      numFailed += testNormalize(OM,S);
    }


    {
      std::vector<int> ind(1); 
      MVT::MvRandom(*S);

      MyOM->stream(Errors) << " projectAndNormalize(): testing on random multivector " << endl;
      numFailed += testProjectAndNormalize(OM,S,X1,X2);
    }


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

      MyOM->stream(Errors) << " projectAndNormalize(): testing [X1 X2]-range multivector against P_X1 P_X2 " << endl;
      numFailed += testProjectAndNormalize(OM,S,X1,X2);
    }


    if (sizeS > 2) {
      MVT::MvRandom(*S);
      RCP<MV> mid = MVT::Clone(*S,1);
      serial_matrix_type c(sizeS,1);
      MVT::MvTimesMatAddMv(ONE,*S,c,ZERO,*mid);
      std::vector<int> ind(1); 
      ind[0] = sizeS-1;
      MVT::SetBlock(*mid,ind,*S);

      MyOM->stream(Errors) << " projectAndNormalize(): testing on rank-deficient multivector " << endl;
      numFailed += testProjectAndNormalize(OM,S,X1,X2);
    }


    if (sizeS > 1) {
      // rank-1
      RCP<MV> one = MVT::Clone(*S,1);
      MVT::MvRandom(*one);
      // put multiple of column 0 in columns 0:sizeS-1
      for (int i=0; i<sizeS; i++) {
        std::vector<int> ind(1); 
        ind[0] = i;
        RCP<MV> Si = MVT::CloneViewNonConst(*S,ind);
        MVT::MvAddMv(SCT::random(),*one,ZERO,*one,*Si);
      }

      MyOM->stream(Errors) << " projectAndNormalize(): testing on rank-1 multivector " << endl;
      numFailed += testProjectAndNormalize(OM,S,X1,X2);
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,cout,success);

  if (numFailed || success==false) {
    if (numFailed) {
      MyOM->stream(Errors) << numFailed << " errors." << endl;
    }
    MyOM->stream(Errors) << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  MyOM->stream(Errors) << "End Result: TEST PASSED" << endl;
  return 0;
}	





static int 
testProjectAndNormalize (RCP< OrthoManager< ST, MV > > OM, 
			 RCP< const MV > S, 
			 RCP< const MV > X1, 
			 RCP< const MV > X2) 
{
  typedef Array< RCP< MV > >::size_type size_type;

  const ST ONE = SCT::one();
  const MT ZERO = SCT::magnitude(SCT::zero());
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

  int numtests;
  numtests = 4;

  // test ortho error before orthonormalizing
  if (X1 != null) {
    MT err = OM->orthogError(*S,*X1);
    sout << "   || <S,X1> || before     : " << err << endl;
  }
  if (X2 != null) {
    MT err = OM->orthogError(*S,*X2);
    sout << "   || <S,X2> || before     : " << err << endl;
  }

  for (int t=0; t<numtests; t++) {

    Array<RCP<const MV> > theX;
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
      // run test
      int ret = OM->projectAndNormalize(*Scopy,theX,C,B);
      sout << "projectAndNormalize() returned rank " << ret << endl;
      if (ret == 0) {
        sout << "   Cannot continue." << endl;
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
        vector<int> ind(ret);
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
        // randomize this data, it should be overwritten
        B->random();
        for (size_type i=0; i<C.size(); i++) {
          C[i]->random();
        }
        // flip the inputs
        theX = tuple( theX[1], theX[0] );
        // run test
        ret = OM->projectAndNormalize(*Scopy,theX,C,B);
        sout << "projectAndNormalize() returned rank " << ret << endl;
        if (ret == 0) {
          sout << "   Cannot continue." << endl;
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
          vector<int> ind(ret);
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
          MT err = OM->orthonormError(*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
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
          MT err = MVDiff(*tmp,*S);
          if (err > ATOL*TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| S_in - X1*C1 - X2*C2 - S_out*B || : " << err << endl;
        }
        // <X1,S> == 0
        if (theX.size() > 0 && theX[0] != null) {
          MT err = OM->orthogError(*theX[0],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <X[0],S> || after      : " << err << endl;
        }
        // <X2,S> == 0
        if (theX.size() > 1 && theX[1] != null) {
          MT err = OM->orthogError(*theX[1],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <X[1],S> || after      : " << err << endl;
        }
      }
    }
    catch (OrthoError e) {
      sout << "   -------------------------------------------         projectAndNormalize() threw exception" << endl;
      sout << "   Error: " << e.what() << endl;
      numerr++;
    }

  } // test for

  MsgType type = Warnings;
  if (numerr>0) type = Errors;
  MyOM->stream(type) << sout.str();
  MyOM->stream(type) << endl;

  return numerr;
}




static int 
testNormalize (RCP< OrthoManager< ST, MV > > OM, 
	       RCP< const MV > S)
{

  const ST ONE = SCT::one();
  const MT ZERO = SCT::magnitude(SCT::zero());
  const int sizeS = MVT::GetNumberVecs(*S);
  int numerr = 0;
  std::ostringstream sout;

  //
  // output tests:
  //   <S_out,S_out> = I
  //   S_in = S_out B
  // 
  // we will loop over an integer specifying the test combinations
  // the bit pattern for the different tests is listed in parenthesis
  //
  // for each of the following, we should test B
  //
  // if hasM:
  // with and without MS  (1)
  //
  // as hasM controls the upper level bits, we need only run test case 0 if hasM==false
  // otherwise, we run test cases 0-1
  //

  int numtests;
  numtests = 1;

  for (int t=0; t<numtests; t++) {

    RCP<serial_matrix_type > B = rcp( new serial_matrix_type(sizeS,sizeS) );

    try {
      // call routine
      // test all outputs for correctness

      // here is where the outputs go
      RCP<MV> Scopy;
      int ret;

      // copies of S,MS
      Scopy = MVT::CloneCopy(*S);
      // randomize this data, it should be overwritten
      B->random();
      // run test
      ret = OM->normalize(*Scopy,B);
      sout << "normalize() returned rank " << ret << endl;
      if (ret == 0) {
        sout << "   Cannot continue." << endl;
        numerr++;
        break;
      }
      // normalize() is only required to return a 
      // basis of rank "ret"
      // this is what we will test:
      //   the first "ret" columns in Scopy
      //   the first "ret" rows in B
      // get pointers to the parts that we want
      if (ret < sizeS) {
        vector<int> ind(ret);
        for (int i=0; i<ret; i++) {
          ind[i] = i;
        }
        Scopy = MVT::CloneViewNonConst(*Scopy,ind);
        B = rcp( new serial_matrix_type(Teuchos::View,*B,ret,sizeS) );
      }

      // test all outputs for correctness
      // S^T M S == I
      {
        MT err = OM->orthonormError(*Scopy);
        if (err > TOL) {
          sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "   || <S,S> - I || after  : " << err << endl;
      }
      // S_in = S_out*B
      {
        RCP<MV> tmp = MVT::Clone(*S,sizeS);
        MVT::MvTimesMatAddMv(ONE,*Scopy,*B,ZERO,*tmp);
        MT err = MVDiff(*tmp,*S);
        if (err > ATOL*TOL) {
          sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| S_in - S_out*B || : " << err << endl;
      }
    }
    catch (OrthoError e) {
      sout << "   -------------------------------------------         normalize() threw exception" << endl;
      sout << "   Error: " << e.what() << endl;
      numerr++;
    }

  } // test for

  MsgType type = Warnings;
  if (numerr>0) type = Errors;
  MyOM->stream(type) << sout.str();
  MyOM->stream(type) << endl;

  return numerr;
}



static int 
testProject (RCP< OrthoManager< ST, MV > > OM, 
	     RCP< const MV > S, 
	     RCP< const MV > X1, 
	     RCP< const MV > X2) 
{
  typedef Array< RCP< MV > >::size_type size_type;

  const ST ONE = SCT::one();
  const int sizeS = MVT::GetNumberVecs(*S);
  const int sizeX1 = MVT::GetNumberVecs(*X1);
  const int sizeX2 = MVT::GetNumberVecs(*X2);
  int numerr = 0;
  std::ostringstream sout;

  //
  // output tests:
  //   <S_out,X1> = 0
  //   <S_out,X2> = 0
  //   S_in = S_out + X1 C1 + X2 C2
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
  // for each of these, we should test 
  // with C1 and C2
  //
  // if hasM:
  // with and without MX1   (1--) 
  // with and without MX2  (1---) 
  // with and without MS  (1----) 
  //
  // as hasM controls the upper level bits, we need only run test cases 0-3 if hasM==false
  // otherwise, we run test cases 0-31
  //

  int numtests;
  numtests = 8;

  // test ortho error before orthonormalizing
  if (X1 != null) {
    MT err = OM->orthogError(*S,*X1);
    sout << "   || <S,X1> || before     : " << err << endl;
  }
  if (X2 != null) {
    MT err = OM->orthogError(*S,*X2);
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
	Array<RCP<MV> > S_outs;
	Array<Array<RCP<serial_matrix_type > > > C_outs;
	RCP<MV> Scopy;

	// copies of S,MS
	Scopy = MVT::CloneCopy(*S);
	// randomize this data, it should be overwritten
	for (size_type i = 0; i < C.size(); ++i) {
	  C[i]->random();
	}
	// run test
	OM->project(*Scopy,theX,C);
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
        // run test
        OM->project(*Scopy,theX,C);
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
          MT err = MVDiff(*tmp,*S);
          if (err > ATOL*TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| S_in - X1*C1 - X2*C2 - S_out || : " << err << endl;
        }
        // <X1,S> == 0
        if (theX.size() > 0 && theX[0] != null) {
          MT err = OM->orthogError(*theX[0],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <X[0],S> || after      : " << err << endl;
        }
        // <X2,S> == 0
        if (theX.size() > 1 && theX[1] != null) {
          MT err = OM->orthogError(*theX[1],*S_outs[o]);
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
          MT err = MVDiff(*S_outs[o1],*S_outs[o2]);
          if (err > TOL) {
            sout << "    vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
        }
      }

    }
    catch (OrthoError e) {
      sout << "   -------------------------------------------         project() threw exception" << endl;
      sout << "   Error: " << e.what() << endl;
      numerr++;
    }

  } // test for

  MsgType type = Warnings;
  if (numerr>0) type = Errors;
  MyOM->stream(type) << sout.str();
  MyOM->stream(type) << endl;

  return numerr;
}



static MT 
MVDiff (const MV& X, 
	const MV& Y) 
{
  const ST ONE = SCT::one();
  const int ncols_X = MVT::GetNumberVecs(X);
  TEST_FOR_EXCEPTION( (MVT::GetNumberVecs(Y) != ncols_X),
		      std::logic_error,
		      "MVDiff: X and Y should have the same number of columns."
		      "  X has " << ncols_X << " column(s) and Y has " 
		      << MVT::GetNumberVecs(Y) << " columns." );
  serial_matrix_type C (ncols_X, ncols_X);

  // tmp := X
  RCP< MV > tmp = MVT::CloneCopy(X);
  // tmp := tmp - Y
  MVT::MvAddMv (-ONE, Y, ONE, *tmp, *tmp);
  // $C := (X - Y)^* \cdot (X - Y)$
  MVT::MvTransMv (ONE, *tmp, *tmp, C);

  // Compute and return $\sum_{j=1}^n \| X(:,j) - Y(:,j) \|_2$, where
  // $n$ is the number of columns in X.
  MT err (0);
  for (int i = 0; i < ncols_X; ++i)
    err += SCT::magnitude (C(i,i));

  return SCT::magnitude (SCT::squareroot (err));
}
