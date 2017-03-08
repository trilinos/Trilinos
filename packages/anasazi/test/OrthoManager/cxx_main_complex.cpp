// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER
//
//  This test is for the OrthoManager interface to ICGSOrthoManager,
//  SVQBOrthoManager and BasicOrthoManager
//
// The matrix used is from MatrixMarket:
// Name: MHD1280B: Alfven Spectra in Magnetohydrodynamics
// Source: Source: A. Booten, M.N. Kooper, H.A. van der Vorst, S. Poedts and J.P. Goedbloed University of Utrecht, the Netherlands
// Discipline: Plasma physics
// URL: http://math.nist.gov/MatrixMarket/data/NEP/mhd/mhd1280b.html
// Size: 1280 x 1280
// NNZ: 22778 entries
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "Trilinos_Util_iohb.h"
#endif

// templated multivector and sparse matrix classes
#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"

using namespace Teuchos;
using namespace Anasazi;
using namespace std;

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;
#else
  typedef double ST;
  // no complex. quit with failure.
int main() {
  std::cout << "Not compiled with complex support." << std::endl;
  std::cout << "End Result: TEST FAILED" << std::endl;
  return -1;
}
#endif

typedef MultiVec<ST>              MV;
typedef Operator<ST>              OP;
typedef MultiVecTraits<ST,MV>    MVT;
typedef OperatorTraits<ST,MV,OP> OPT;
typedef ScalarTraits<ST>         SCT;
typedef SCT::magnitudeType           MT;

// this is the tolerance that all tests are performed against
const MT TOL = 1.0e-12;
const MT ATOL = 10;

// declare an output manager for handling local output
RCP< Anasazi::BasicOutputManager<ST> > MyOM;

// some forward declarations
int testProject(RCP<OrthoManager<ST,MV> > OM, RCP<const MV> S, RCP<const MV> X1, RCP<const MV> X2);
int testNormalize(RCP<OrthoManager<ST,MV> > OM, RCP<const MV> S);
int testProjectAndNormalize(RCP<OrthoManager<ST,MV> > OM, RCP<const MV> S, RCP<const MV> X1, RCP<const MV> X2);

MT MVDiff(const MV &X, const MV &Y);

int main(int argc, char *argv[]) 
{
  
  const ST ONE = SCT::one();
  const MT ZERO = SCT::magnitude(SCT::zero());
  bool verbose = false;
  int numFailed = 0;
  bool debug = false;
  string ortho = "SVQB";
  std::string filename;
  int dim = 100;
  int sizeS  = 5;
  int sizeX1 = 11; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  int sizeX2 = 13; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1

  bool success = true;
  try {

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("ortho",&ortho,"Which ortho manager: SVQB or Basic or ICGS");
    cmdp.setOption("dim",&dim,"Controls the size of multivectors.");
    cmdp.setOption("sizeS",&sizeS,"Controls the width of the input multivector.");
    cmdp.setOption("sizeX1",&sizeX1,"Controls the width of the first basis.");
    cmdp.setOption("sizeX2",&sizeX2,"Controls the width of the second basis.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }

    // below we will assume that sizeX1 > sizeX2
    // this does not effect our testing, since we will test P_{X1,Y1} P_{X2,Y2} as well as P_{X2,Y2} P_{X1,Y1}
    // however, is does allow us to simplify some logic
    if (sizeX1 < sizeX2) {
      swap(sizeX1,sizeX2);
    }

    // instantiate the output manager
    MyOM = rcp( new BasicOutputManager<ST>() );
    if (verbose) {
      // output in this driver will be sent to Anasazi::Warnings
      MyOM->setVerbosity(Anasazi::Warnings);
    }

    // Output Anasazi version
    MyOM->stream(Anasazi::Warnings) << Anasazi_Version() << endl << endl;

    //  Problem information
    RCP< const MyBetterOperator<ST> > M;
    if (filename != "") {
      int dim2,nnz;
      double *dvals;
      int *colptr,*rowind;
      nnz = -1;
      int info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,
          &colptr,&rowind,&dvals);
      TEUCHOS_TEST_FOR_EXCEPTION(info == 0 || nnz < 0, 
          std::runtime_error,
          "Error reading file " << filename << ". info == " << info << " and nnz == " << nnz);
      // Convert interleaved doubles to complex values
      std::vector<ST> cvals(nnz);
      for (int ii=0; ii<nnz; ii++) {
        cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
      }
      // Build the problem matrix
      M = rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,&cvals[0]) );
    }

    // Create ortho managers
    RCP<OrthoManager<ST,MV> > OM;
    if (ortho == "SVQB") {
      OM = rcp( new SVQBOrthoManager<ST,MV,OP>(M) );
    }
    else if (ortho == "Basic") {
      OM = rcp( new BasicOrthoManager<ST,MV,OP>(M) );
    }
    else if (ortho == "ICGS") {
      OM = rcp( new ICGSOrthoManager<ST,MV,OP>(M) );
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,invalid_argument,"Command line parameter \"ortho\" must be \"SVQB\" or \"Basic\".");
    }


    // multivector to spawn off of
    RCP<MV> S = rcp( new MyMultiVec<ST>(dim, sizeS) );

    // create X1, X2
    // they must be M-orthonormal and mutually M-orthogonal
    MyOM->stream(Errors) << " Generating X1,X2 for testing... " << endl;
    RCP<MV> X1  = MVT::Clone(*S,sizeX1),
      X2  = MVT::Clone(*S,sizeX2);
    {
      int dummy;
      MT err;
      // X1
      MVT::MvRandom(*X1);
      dummy = OM->normalize(*X1);
      TEUCHOS_TEST_FOR_EXCEPTION(dummy != sizeX1, std::runtime_error, 
          "normalize(X1) returned rank " << dummy << " from " 
          << sizeX1 << " vectors. Cannot continue.");
      err = OM->orthonormError(*X1);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "normalize(X1) did meet tolerance: orthonormError(X1) == " << err);
      MyOM->stream(Warnings) << "   || <X1,X1> - I || : " << err << endl;
      // X2
      MVT::MvRandom(*X2);
      dummy = OM->projectAndNormalize(*X2,tuple<RCP<const MV> >(X1));
      TEUCHOS_TEST_FOR_EXCEPTION(dummy != sizeX2, std::runtime_error, 
          "projectAndNormalize(X2,X1) returned rank " << dummy << " from " 
          << sizeX2 << " vectors. Cannot continue.");
      err = OM->orthonormError(*X2);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "projectAndNormalize(X2,X1) did not meet tolerance: orthonormError(X2) == " << err);
      err = OM->orthogError(*X2,*X1);
      MyOM->stream(Warnings) << "   || <X2,X2> - I || : " << err << endl;
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "projectAndNormalize(X2,X1) did not meet tolerance: orthogError(X2,X1) == " << err);
      MyOM->stream(Warnings) << "   || <X2,X1> ||     : " << err << endl;
    }
    MyOM->stream(Warnings) << endl;


    {
      // just a random multivector
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
      SerialDenseMatrix<int,ST> C1(sizeX1,sizeS), C2(sizeX2,sizeS);
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
      SerialDenseMatrix<int,ST> c(sizeS,1);
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
      SerialDenseMatrix<int,ST> C1(sizeX1,sizeS), C2(sizeX2,sizeS);
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
      SerialDenseMatrix<int,ST> c(sizeS,1);
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





////////////////////////////////////////////////////////////////////////////
int testProjectAndNormalize(RCP<OrthoManager<ST,MV> > OM, 
                            RCP<const MV> S, 
                            RCP<const MV> X1, RCP<const MV> X2) {

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
    RCP<SerialDenseMatrix<int,ST> > B = rcp( new SerialDenseMatrix<int,ST>(sizeS,sizeS) );
    Array<RCP<SerialDenseMatrix<int,ST> > > C;
    if ( (t && 3) == 0 ) {
      // neither <X1,Y1> nor <X2,Y2>
      // C, theX and theY are already empty
    }
    else if ( (t && 3) == 1 ) {
      // X1
      theX = tuple(X1);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)) );
    }
    else if ( (t && 3) == 2 ) {
      // X2
      theX = tuple(X2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }
    else {
      // X1 and X2, and the reverse.
      theX = tuple(X1,X2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)), 
                 rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }

    try {
      // call routine
      // if (t && 3) == 3, {
      //    call with reversed input: X2 X1
      // }
      // test all outputs for correctness
      // test all outputs for equivalence

      // here is where the outputs go
      Array<RCP<const MV> > S_outs;
      Array<Array<RCP<SerialDenseMatrix<int,ST> > > > C_outs;
      Array<RCP<SerialDenseMatrix<int,ST> > > B_outs;
      RCP<MV> Scopy;
      Array<int> ret_out;

      // copies of S,MS
      Scopy = MVT::CloneCopy(*S);
      // randomize this data, it should be overwritten
      B->random();
      for (unsigned int i=0; i<C.size(); i++) {
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
        S_outs.push_back( MVT::CloneView(*Scopy,ind) );
        B_outs.push_back( rcp( new SerialDenseMatrix<int,ST>(Teuchos::Copy,*B,ret,sizeS) ) );
      }
      else {
        S_outs.push_back( Scopy );
        B_outs.push_back( rcp( new SerialDenseMatrix<int,ST>(*B) ) );
      }
      C_outs.push_back( Array<RCP<SerialDenseMatrix<int,ST> > >(0) );
      if (C.size() > 0) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[0]) ) );
      }
      if (C.size() > 1) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[1]) ) );
      }

      // do we run the reversed input?
      if ( (t && 3) == 3 ) {
        // copies of S,MS
        Scopy = MVT::CloneCopy(*S);
        // randomize this data, it should be overwritten
        B->random();
        for (unsigned int i=0; i<C.size(); i++) {
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
          S_outs.push_back( MVT::CloneView(*Scopy,ind) );
          B_outs.push_back( rcp( new SerialDenseMatrix<int,ST>(Teuchos::Copy,*B,ret,sizeS) ) );
        }
        else {
          S_outs.push_back( Scopy );
          B_outs.push_back( rcp( new SerialDenseMatrix<int,ST>(*B) ) );
        }
        C_outs.push_back( Array<RCP<SerialDenseMatrix<int,ST> > >() );
        // reverse the Cs to compensate for the reverse projectors
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[1]) ) );
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[0]) ) );
        // flip the inputs back
        theX = tuple( theX[1], theX[0] );
      }


      // test all outputs for correctness
      for (unsigned int o=0; o<S_outs.size(); o++) {
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
    catch (const OrthoError &e) {
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



////////////////////////////////////////////////////////////////////////////
int testNormalize(RCP<OrthoManager<ST,MV> > OM, RCP<const MV> S)
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

    // allocate the data separately from B, to allow it to persist as B is delete/reallocated
    Teuchos::ArrayRCP<ST> Bdata = Teuchos::arcp<ST>(sizeS*sizeS);
    RCP<SerialDenseMatrix<int,ST> > B = rcp( new SerialDenseMatrix<int,ST>(Teuchos::View,Bdata.getRawPtr(),sizeS,sizeS,sizeS) );

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
        B = rcp( new SerialDenseMatrix<int,ST>(Teuchos::View,Bdata.getRawPtr(),ret,ret,sizeS) );
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
    catch (const OrthoError &e) {
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



////////////////////////////////////////////////////////////////////////////
int testProject(RCP<OrthoManager<ST,MV> > OM, 
                   RCP<const MV> S, 
                   RCP<const MV> X1, RCP<const MV> X2) {

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

  for (int t=0; t<numtests; t++) {

    Array<RCP<const MV> > theX;
    Array<RCP<SerialDenseMatrix<int,ST> > > C;
    if ( (t && 3) == 0 ) {
      // neither X1 nor X2
      // C and theX are already empty
    }
    else if ( (t && 3) == 1 ) {
      // X1
      theX = tuple(X1);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)) );
    }
    else if ( (t && 3) == 2 ) {
      // X2
      theX = tuple(X2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }
    else {
      // X1 and X2, and the reverse.
      theX = tuple(X1,X2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)), 
                 rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
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
      Array<Array<RCP<SerialDenseMatrix<int,ST> > > > C_outs;
      RCP<MV> Scopy;

      // copies of S,MS
      Scopy = MVT::CloneCopy(*S);
      // randomize this data, it should be overwritten
      for (unsigned int i=0; i<C.size(); i++) {
        C[i]->random();
      }
      // run test
      OM->project(*Scopy,theX,C);
      // we allocate S and MS for each test, so we can save these as views
      // however, save copies of the C
      S_outs.push_back( Scopy );
      C_outs.push_back( Array<RCP<SerialDenseMatrix<int,ST> > >(0) );
      if (C.size() > 0) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[0]) ) );
      }
      if (C.size() > 1) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[1]) ) );
      }

      // do we run the reversed input?
      if ( (t && 3) == 3 ) {
        // copies of S,MS
        Scopy = MVT::CloneCopy(*S);
        // randomize this data, it should be overwritten
        for (unsigned int i=0; i<C.size(); i++) {
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
        C_outs.push_back( Array<RCP<SerialDenseMatrix<int,ST> > >() );
        // reverse the Cs to compensate for the reverse projectors
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[1]) ) );
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[0]) ) );
        // flip the inputs back
        theX = tuple( theX[1], theX[0] );
      }

      // test all outputs for correctness
      for (unsigned int o=0; o<S_outs.size(); o++) {
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
      for (unsigned int o1=0; o1<S_outs.size(); o1++) {
        for (unsigned int o2=o1+1; o2<S_outs.size(); o2++) {
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
    catch (const OrthoError &e) {
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



MT MVDiff(const MV &X, const MV &Y) {
  const ST ONE = SCT::one();
  const int sizeX = MVT::GetNumberVecs(X);
  SerialDenseMatrix<int,ST> xTmx(sizeX,sizeX);

  // tmp <- X
  RCP<MV> tmp = MVT::CloneCopy(X);

  // tmp <- tmp - Y
  MVT::MvAddMv(-ONE,Y,ONE,*tmp,*tmp);

  MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
  MT err = 0;
  for (int i=0; i<sizeX; i++) {
    err += SCT::magnitude(xTmx(i,i));
  }
  return SCT::magnitude(SCT::squareroot(err));
}
