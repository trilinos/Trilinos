// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
//
//  This test is for the SVQBOrthoManager
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziSolverUtils.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

using namespace Teuchos;
using namespace Anasazi;
using namespace std;

typedef double                       ST;
typedef Epetra_MultiVector           MV;
typedef Epetra_Operator              OP;
typedef MultiVecTraits<double,MV>    MVT;
typedef OperatorTraits<double,MV,OP> OPT;
typedef ScalarTraits<double>         SCT;
typedef SCT::magnitudeType           MT;

const ST ONE = SCT::one();
const MT ZERO = SCT::magnitude(SCT::zero());

// this is the tolerance that all tests are performed against
const MT TOL = 1.0e-13;
const MT ATOL = 100;

// declare an output manager for handling local output
RCP< Anasazi::BasicOutputManager<ST> > MyOM;

// some forward declarations
int testProjectAndNormalizeGen(RCP<GenOrthoManager<ST,MV,OP> > OM, RCP<const MV> S, RCP<const MV> X1, RCP<const MV> Y1, RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho);
int testProjectGen(RCP<GenOrthoManager<ST,MV,OP> > OM, RCP<const MV> S, RCP<const MV> X1, RCP<const MV> Y1, RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho);

MT MVDiff(const MV &X, const MV &Y);

int main(int argc, char *argv[]) 
{
  
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = false;
  int numFailed = 0;
  bool debug = false;
  string ortho = "SVQB";
  bool useMass = true;
  int numElements = 100;
  int sizeS  = 5;
  int sizeX1 = 11; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  int sizeX2 = 13; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("ortho",&ortho,"Which ortho manager: SVQB or Basic");
  cmdp.setOption("numElements",&numElements,"Controls the size of the mass matrix.");
  cmdp.setOption("useMass","noMass",&useMass,"Use a mass matrix for inner product.");
  cmdp.setOption("sizeS",&sizeS,"Controls the width of the input multivector.");
  cmdp.setOption("sizeX1",&sizeX1,"Controls the width of the first basis.");
  cmdp.setOption("sizeX2",&sizeX2,"Controls the width of the second basis.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
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

  // Problem information
  RCP<const Epetra_CrsMatrix> M;
  if (useMass) {
    const int space_dim = 1;
    std::vector<ST> brick_dim( space_dim );
    brick_dim[0] = 1.0;
    std::vector<int> elements( space_dim );
    elements[0] = numElements;

    // Create problem
    RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
    // Get the mass matrix
    RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  }

  // Create ortho managers
  RCP<GenOrthoManager<ST,MV,OP> > OM;
  if (ortho == "SVQB") {
    OM = rcp( new SVQBOrthoManager<ST,MV,OP>(M) );
  }
  else if (ortho == "Basic") {
    OM = rcp( new BasicOrthoManager<ST,MV,OP>(M) );
  }
  else {
    TEST_FOR_EXCEPTION(true,invalid_argument,"Command line parameter \"ortho\" must be \"SVQB\" or \"Basic\".");
  }

  // multivector to spawn off of
  RCP<MV> S = rcp( new Epetra_MultiVector(M->OperatorDomainMap(), sizeS) );

  // create X1, Y1, X2, Y2
  // need <X1,Y2> = 0 = <X2,Y1>... we can use the ortho manager to give us this
  // want a set so that <X1b,Y1> = I = <X2b,Y2>... do this via modifications of X1,X2
  // also want <Y1,Y1> = I = <Y2,Y2> so that we can test P_Y1,Y1 and P_Y2,Y2
  MyOM->stream(Errors) << " Generating Y1,Y2 for project() : testing... " << endl;
  RCP<MV> X1  = MVT::Clone(*S,sizeX1),
          Y1  = MVT::Clone(*S,sizeX1),
          X1b = MVT::Clone(*S,sizeX1),
          X2  = MVT::Clone(*S,sizeX2),
          Y2  = MVT::Clone(*S,sizeX2),
          X2b = MVT::Clone(*S,sizeX2);
  {
    int dummy;
    MT err;
    SerialDenseMatrix<int,ST> yTx1(sizeX1,sizeX1), yTx2(sizeX2,sizeX2);
    Array<int> yTx1_piv(sizeX1), yTx2_piv(sizeX2);
    Teuchos::LAPACK<int,ST> lapack;
    TEST_FOR_EXCEPTION(sizeX1 < sizeX2,std::logic_error,"Internal logic error: sizeX1 < sizeX2.");
    Array<ST> LUwork(sizeX1);
    // Y1
    MVT::MvRandom(*Y1);
    dummy = OM->normalize(*Y1);
    TEST_FOR_EXCEPTION(dummy != sizeX1, std::runtime_error, 
        "normalize(Y1) returned rank " << dummy << " from " 
        << sizeX1 << " vectors. Cannot continue.");
    err = OM->orthonormError(*Y1);
    TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
        "normalize(Y1) did meet tolerance: orthonormError(Y1) == " << err);
    MyOM->stream(Warnings) << "   || <Y1,Y1> - I || : " << err << endl << endl;
    // Y2
    MVT::MvRandom(*Y2);
    dummy = OM->normalize(*Y2);
    TEST_FOR_EXCEPTION(dummy != sizeX1, std::runtime_error, 
        "normalize(Y2) returned rank " << dummy << " from " 
        << sizeX1 << " vectors. Cannot continue.");
    err = OM->orthonormError(*Y2);
    TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
        "normalize(Y2) did not meet tolerance: orthonormError(Y2) == " << err);
    MyOM->stream(Warnings) << "   || <Y2,Y2> - I || : " << err << endl << endl;
    // X2 ortho to Y1
    MVT::MvRandom(*X2);
    OM->project(*X2,*Y1);
    err = OM->orthogError(*X2,*Y1);
    TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
        "project(X2,Y1) did not meet tolerance: orthog(X2,Y1) == " << err);
    MyOM->stream(Warnings) << "   || <X2,Y1> ||     : " << err << endl << endl;
    // X1 ortho to Y2
    MVT::MvRandom(*X1);
    OM->project(*X1,*Y2);
    err = OM->orthogError(*X1,*Y2);
    TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
        "project(X1,Y2) did not meet tolerance: orthog(X1,Y2) == " << err);
    MyOM->stream(Warnings) << "   || <X1,Y2> ||     : " << err << endl << endl;
    // Compute X1b so that <X1b,Y1> = I
    // Compute LU factorization of <Y1,X1> and use it to compute explicit inverse of <Y1,X1>
    OM->innerProd(*Y1,*X1,yTx1);
    lapack.GETRF(yTx1.numRows(),yTx1.numCols(),yTx1.values(),yTx1.stride(),&yTx1_piv.front(),&dummy);
    TEST_FOR_EXCEPTION(dummy != 0,std::logic_error, "Error computing LU factorization of <Y1,X1>.");
    lapack.GETRI(yTx1.numRows(),yTx1.values(),yTx1.stride(),&yTx1_piv.front(),&LUwork.front(),LUwork.size(),&dummy);
    TEST_FOR_EXCEPTION(dummy != 0,std::logic_error, "Error computing LU inverse of <Y1,X1>.");
    // Set X1b=X1*inv(<Y1,X1>), so that <Y1,X1b> = <Y1,X1*inv(<Y1,X1>)> = <Y1,X1>*inv(<Y1,X1>) = I
    MVT::MvTimesMatAddMv(ONE,*X1,yTx1,ZERO,*X1b);
    // Test that this was successful
    OM->innerProd(*Y1,*X1b,yTx1);
    for (int i=0; i<sizeX1; i++) {
      yTx1(i,i) -= ONE;
    }
    err = yTx1.normFrobenius();
    TEST_FOR_EXCEPTION(err > TOL,std::runtime_error, "Failed to make <X1b,Y1> = I.");
    // Compute X2b so that <X2b,Y2> = I
    // Compute LU factorization of <Y2,X2> and use it to compute explicit inverse of <Y2,X2>
    OM->innerProd(*Y2,*X2,yTx2);
    lapack.GETRF(yTx2.numRows(),yTx2.numCols(),yTx2.values(),yTx2.stride(),&yTx2_piv.front(),&dummy);
    TEST_FOR_EXCEPTION(dummy != 0,std::logic_error, "Error computing LU factorization of <Y2,X2>.");
    lapack.GETRI(yTx2.numRows(),yTx2.values(),yTx2.stride(),&yTx2_piv.front(),&LUwork.front(),LUwork.size(),&dummy);
    TEST_FOR_EXCEPTION(dummy != 0,std::logic_error, "Error computing LU inverse of <Y2,X2>.");
    // Set X2b=X2*inv(<Y2,X2>), so that <Y2,X2b> = <Y2,X2*inv(<Y2,X2>)> = <Y2,X2>*inv(<Y2,X2>) = I
    MVT::MvTimesMatAddMv(ONE,*X2,yTx2,ZERO,*X2b);
    // Test that this was successful
    OM->innerProd(*Y2,*X2b,yTx2);
    for (int i=0; i<sizeX2; i++) {
      yTx2(i,i) -= ONE;
    }
    err = yTx2.normFrobenius();
    TEST_FOR_EXCEPTION(err > TOL,std::runtime_error, "Failed to make <X2b,Y2> = I.");
  }

  {
    MVT::MvRandom(*S);
    
    MyOM->stream(Errors) << " projectGen(): testing on random multivector " << endl;
    numFailed += testProjectGen(OM,S,X1,Y1,X2,Y2,false);
    MyOM->stream(Errors) << " projectGen(biOrtho): testing on random multivector " << endl;
    numFailed += testProjectGen(OM,S,X1b,Y1,X2b,Y2,true);
  }

  {
    std::vector<int> ind(1); 
    MVT::MvRandom(*S);
    
    MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on random multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
    MyOM->stream(Errors) << " projectAndNormalizeGen(biOrtho): testing on random multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
  }

  {
    // finish: what to do here...
    /*
    SerialDenseMatrix<int,ST> B(sizeQ,sizeS);
    B.random();
    MVT::MvTimesMatAddMv(ONE,*Q1,B,ZERO,*S);
    
    MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on Q-range multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
    MyOM->stream(Errors) << " projectAndNormalizeGen(biOrtho): testing on Q-range multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
    */
  }

  {
    std::vector<int> ind(1); 
    MVT::MvRandom(*S);
    ind[0] = 1;
    RCP<MV> mid = MVT::CloneCopy(*S,ind);
    ind[0] = 2;
    MVT::SetBlock(*mid,ind,*S);
    
    MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on rank-deficient multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
    MyOM->stream(Errors) << " projectAndNormalizeGen(biOrtho): testing on rank-deficient multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
  }

  {
    MVT::MvRandom(*S);
    std::vector<int> ind(1); 
    // get column 0
    ind[0] = 0;
    RCP<MV> mid = MVT::CloneCopy(*S,ind);
    // put column 0 in columns 1:sizeS-1
    for (int i=1; i<sizeS; i++) {
      ind[0] = i;
      MVT::SetBlock(*mid,ind,*S);
    }
    
    MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on rank-1 multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
    MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on rank-1 multivector " << endl;
    numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (numFailed) {
    MyOM->stream(Errors) << numFailed << " errors." << endl;
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
int testProjectAndNormalizeGen(RCP<GenOrthoManager<ST,MV,OP> > OM, 
                               RCP<MV> X, 
                               RCP<const MV> X1, RCP<const MV> Y1,
                               RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho) {

  const int sizeX = MVT::GetNumberVecs(*X);
  const int sizeX1 = MVT::GetNumberVecs(*X1);
  const int sizeX2 = MVT::GetNumberVecs(*X2);
  MT err;
  int rank;
  RCP<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  int numerr = 0;
  bool hasM = (OM->getOp() != null);
  std::ostringstream sout;
  RCP<MV> xcopy, mxcopy, smlOX;
  RCP<SerialDenseMatrix<int,ST> > C, smlC, R, smlR, newR;
  bool warning = false;

  //
  // output tests:
  //   <S_out,S_out> = I
  //   <S_out,Y1> = 0
  //   <S_out,Y2> = 0
  //   S_in = S_out B + X1 C1 + X2 C2
  // 
  // we will loop over an integer specifying the test combinations
  // the bit pattern for the different tests is listed in parenthesis
  //
  // for the projectors, test the following combinations:
  // none                        (00)
  // P_{X1,Y1}                   (01)
  // P_{X2,Y2}                   (10)
  // P_{X1,Y1} P_{X2,Y2}         (11)
  // P_{X2,Y2} P_{X1,Y1}         (11)
  // the latter two should be tested to give the same answer
  // this relies on <X1,Y2> = 0 = <X2,Y1>, ensured in the main() routine
  //
  // for each of these, we should test 
  // with C1, C2 and B
  // with and without isBiortho (if isBiortho==true)  (1--)
  //
  // comparison should be made against projectAndNormalizeMat() if X1==Y1 and X2==Y2 and isBiortho==true
  //
  // if hasM:
  // with and without MX1,MY1   (1---) 
  // with and without MX2,MY2  (1----) 
  // with and without MS      (1-----) 
  //
  // as hasM controls the upper level bits, we need only run test cases 0-7 if hasM==false
  // otherwise, we run test cases 0-63
  //

  bool hasM = (OM->getOp() != Teuchos::null);
  int numtests;
  RCP<const MV> MX1, MX2, MY1, MY2; 
  RCP<MV> MS;
  if (hasM) {
    MX1 = MVT::Clone(*S,sizeX1);
    MY1 = MVT::Clone(*S,sizeX1);
    MX2 = MVT::Clone(*S,sizeX2);
    MY2 = MVT::Clone(*S,sizeX2);
    MS  = MVT::Clone(*S,sizeS);
    OPT::Apply(*(OM->getOp()),*X1,const_cast<MV>(*MX1));
    OPT::Apply(*(OM->getOp()),*Y1,const_cast<MV>(*MY1));
    OPT::Apply(*(OM->getOp()),*X2,const_cast<MV>(*MX2));
    OPT::Apply(*(OM->getOp()),*Y2,const_cast<MV>(*MY2));
    OPT::Apply(*(OM->getOp()),*S ,*MS);
    numtests = 64;
  }
  else {
    numtests = 7;
  }

  /*
  // test ortho error before orthonormalizing
  err = OM->orthonormError(*X);
  sout << "   || X^T M X - I ||_F before : " << err << endl;
  err = OM->orthogError(*X,*Q);
  sout << "   || Q^T M X ||_F before     : " << err << endl;
  */

  for (int t=0; t<numtests; t++) {

    // make a copy of S to work from
    Scopy = MVT::CloneCopy(*S);

    // pointers to simplify calls below
    RCP<const MV> lclMX1;
    RCP<const MV> lclMY1;
    RCP<const MV> lclMX2;
    RCP<const MV> lclMY2;
    RCP<MV> lclMS;
    if ( t & 8 ) {
      lclMX1 = MX1;
    }
    if ( t & 16 ) {
      lclMY1 = MY1;
    }
    if ( t & 32 ) {
      lclMX2 = MX2;
    }
    if ( t & 64 ) {
      lclMY2 = MY2;
    }
    if ( t & 128 ) {
      lclMS = MS;
    }

    Array<RCP<const MV> > theX, theY;
    RCP<SerialDenseMatrix<int,ST> > B = rcp( new SerialDenseMatrix<int,ST>(sizeS,sizeS) );
    R->random();
    Array<RCP<SerialDenseMatrix<int,ST> > > C;
    if ( (t && 3) == 0 ) {
      // neither <X1,Y1> nor <X2,Y2>
      // C, theX and theY are already empty
    }
    else if ( (t && 3) == 1 ) {
      // <X1,Y1>
      theX = tuple(X1);
      theY = tuple(Y1);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)) );
    }
    else if ( (t && 3) == 2 ) {
      // <X2,Y2>
      theX = tuple(X2);
      theY = tuple(Y2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }
    else {
      // <X1,Y2> and <X2,Y2>, and the reverse.
      theX = tuple(X1,X2);
      theY = tuple(Y1,Y2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)), 
                 rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }

    bool localIsBiortho = isBiortho;
    if (isBiortho) {
    }
    else {
    }

    try {
      // call routine
      // test output
      // if (t && 3) == 3, {
      //    call reversed routine, test output, and compare outputs from previous
      // }

      rank = OM->projectAndNormalizeGen(
          *xcopy,theX,theY,
          tuple< RCP< SerialDenseMatrix<int,ST> > >(C),R, mxcopy
          );

      sout << "projectAndNormalize() returned rank " << rank << endl;

      ind.resize(rank);
      for (int i=0; i<rank; i++) {
        ind[i] = i;
      }
      smlX = MVT::CloneView(*xcopy,ind); 
      smlOX = MVT::CloneView(*X,ind);
      if (mxcopy != null) {
        smlMX = MVT::CloneView(*mxcopy,ind); 
      }
      else {
        smlMX = null;
      }

      // MX == M*X
      if (smlMX != null) {
        tmp = MVT::Clone(*smlX,rank);
        OPT::Apply(*(OM->getOp()),*smlX,*tmp);
        err = MVDiff(*tmp,*smlMX);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| MX - M*X ||_F           : " << err << endl;
      }
      // X = xcopy*R
      if (R != null) {
        smlR = rcp( new SerialDenseMatrix<int,ST>(Teuchos::Copy,*R,rank,rank) );
        newR = rcp( new SerialDenseMatrix<int,ST>(rank,rank) );
        OM->innerProd(*smlX,*smlOX,*newR);
        *newR -= *smlR;
        err = newR->normFrobenius();
        if (err > ATOL*TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        }
        sout << "  " << t << "|| newX'*M*X-R ||_F          : " << err << endl;
      }
      // X^T M X == I
      err = OM->orthonormError(*smlX);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "   || X^T M X - I ||_F after  : " << err << endl;
      // X = Q*C + xcopy*R
      if (C != null) {
        smlR = rcp( new SerialDenseMatrix<int,ST>(Teuchos::Copy,*R,rank,rank) );
        smlC = rcp( new SerialDenseMatrix<int,ST>(Teuchos::Copy,*C,sizeQ,rank) );
        tmp = MVT::Clone(*smlX,rank);
        MVT::MvTimesMatAddMv(ONE,*smlX,*smlR,ZERO,*tmp);
        MVT::MvTimesMatAddMv(ONE,*Q,*smlC,ONE,*tmp);
        err = MVDiff(*tmp,*smlOX);
        if (err > ATOL*TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        }
        sout << "  " << t << "|| X-Q*C-newX*R ||_F       : " << err << endl;
      }
      // Q^T M X == I
      err = OM->orthogError(*smlX,*Q);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| Q^T M X ||_F after      : " << err << endl;

    }
    catch (OrthoError e) {
      sout << "   -------------------------------------------         projectAndNormalize() Failed" << endl;
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



//////////////////////////////////////////////////////////////////////
int testProjectGen(RCP<GenOrthoManager<ST,MV,OP> > OM, 
                RCP<const MV> X, RCP<const MV> Q1, RCP<const MV> Q2) {

  const int sizeX = MVT::GetNumberVecs(*X);
  const int sizeQ = MVT::GetNumberVecs(*Q);
  MT err;
  RCP<MV> tmp;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;
  std::ostringstream sout;
  RCP<MV> xcopy, mxcopy;
  RCP<SerialDenseMatrix<int,ST> > C;
  bool warning = false;

  int numtests;
  RCP<MV> MX;
  if (hasM) {
    numtests = 4;
    MX = MVT::Clone(*X,sizeX);
    OPT::Apply(*(OM->getOp()),*X,*MX);
  }
  else {
    MX = MVT::CloneCopy(*X);
    numtests = 2;
  }

  // test ortho error before orthonormalizing
  err = OM->orthogError(*Q,*X);
  sout << "   || Q^T M X ||_F before     : " << err << endl;

  for (int t=0; t<numtests; t++) {

    xcopy = MVT::CloneCopy(*X);

    if (t & 2) {
      mxcopy = MVT::CloneCopy(*MX);
    }
    else {
      mxcopy = null;
    }

    if (t & 1) {
      C = rcp( new SerialDenseMatrix<int,ST>(sizeQ,sizeX) );
      // whatever is here should be overwritten
      C->random();
    }
    else {
      C = null;
    }

    try {
      OM->projectGen(*xcopy, tuple<RCP<const MV> >(Q), 
                  tuple<RCP<SerialDenseMatrix<int,ST> > >(C), mxcopy
                  );
      // MX == M*X
      if (mxcopy != null) {
        tmp = MVT::CloneCopy(*xcopy);
        OPT::Apply(*(OM->getOp()),*xcopy,*tmp);
        err = MVDiff(*tmp,*mxcopy);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| MX - M*X ||_F           : " << err << endl;
      }
      // X = Q*C + xcopy
      if (C != null) {
        tmp = MVT::CloneCopy(*xcopy);
        MVT::MvTimesMatAddMv(ONE,*Q,*C,ONE,*tmp);
        err = MVDiff(*tmp,*X);
        if (err > ATOL*TOL) {
          warning = true;
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! warning!" << endl;
        }
        sout << "  " << t << "|| X-Q*C-newX ||_F         : " << err << endl;
      }
      // Q^T M X == I
      err = OM->orthogError(*xcopy,*Q);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| Q^T M X ||_F after      : " << err << endl;
    }
    catch (OrthoError e) {
      sout << "   -------------------------------------------         project() Failed" << endl;
      sout << "   Error: " << e.what() << endl;
      numerr++;
    }

  } // for test

  MsgType type = Warnings;
  if (numerr>0 || warning) type = Errors;
  MyOM->stream(type) << sout.str();
  MyOM->stream(type) << endl;

  return numerr;
}




MT MVDiff(const MV &X, const MV &Y) {
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
