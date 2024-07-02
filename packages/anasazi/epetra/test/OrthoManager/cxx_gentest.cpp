// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the GenOrthoManager interface to ICGSOrthoManager
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziSolverUtils.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"
#include "MySDMHelpers.hpp"

using namespace Teuchos;
using namespace Anasazi;
using std::cout;
using std::endl;
using std::swap;

typedef double                       ST;
typedef Epetra_MultiVector           MV;
typedef Epetra_Operator              OP;

// this is the tolerance that all tests are performed against
const double TOL = 1.0e-10;
const double ATOL = 100;

// declare an output manager for handling local output
RCP< Anasazi::BasicOutputManager<ST> > MyOM;

// some forward declarations
int testProjectAndNormalizeGen(RCP<GenOrthoManager<ST,MV,OP> > OM, RCP<const MV> S, RCP<const MV> X1, RCP<const MV> Y1, RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho);
int testProjectGen(RCP<GenOrthoManager<ST,MV,OP> > OM, RCP<const MV> S, RCP<const MV> X1, RCP<const MV> Y1, RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho);

double MVDiff(const MV &X, const MV &Y);

int main(int argc, char *argv[])
{
  typedef MultiVecTraits<double,MV>    MVT;
  //typedef OperatorTraits<double,MV,OP> OPT;
  typedef ScalarTraits<double>         SCT;
  typedef SCT::magnitudeType           MT;

  const ST ONE = SCT::one();
  const MT ZERO = SCT::magnitude(SCT::zero());
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
  bool useMass = true;
  int dim = 100;
  int sizeS  = 3;
  int sizeX1 = 6; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1
  int sizeX2 = 7; // MUST: sizeS + sizeX1 + sizeX2 <= elements[0]-1

  bool success = true;
  try {

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("dim",&dim,"Controls the size of multivectors.");
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
    RCP<ModalProblem> testCase;
    if (useMass) {
      const int space_dim = 1;
      std::vector<ST> brick_dim( space_dim );
      brick_dim[0] = 1.0;
      std::vector<int> elements( space_dim );
      elements[0] = dim+1;

      // Create problem
      testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
      // Get the mass matrix
      M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
    }

    // Create ortho managers
    RCP<GenOrthoManager<ST,MV,OP> > OM = rcp( new ICGSOrthoManager<ST,MV,OP>(M) );

    // multivector to spawn off of
    RCP<MV> S;
    if (useMass) {
      S = rcp( new Epetra_MultiVector(M->OperatorDomainMap(),sizeS) );
    }
    else {
      Epetra_Map map(dim,0,Comm);
      S = rcp( new Epetra_MultiVector(map, sizeS) );
    }

    // create X1, Y1, X2, Y2
    //
    // need <X1,Y2> = 0 = <X2,Y1>... we can use the ortho manager to give us this
    //
    // want a set so that <X1b,Y1> = I = <X2b,Y2>... do this via modifications of X1,X2
    //
    // need <Y1,X1> s.p.d., so let X1 = M*Y1, so that
    // <Y1,X1> = Y1'*M*M*Y1, which is a non-trivial s.p.d. matrix
    // same for <X2,Y2>
    //
    // Y1 cannot be orthonormal, or else <Y1,X1b>==I implies that X1b==Y1
    // same with Y2,X2b
    // therefore, we cannot call
    MyOM->stream(Errors) << " Generating Y1,Y2 for project() : testing... " << endl;
    RCP<MV> X1  = MVT::Clone(*S,sizeX1),
      Y1  = MVT::Clone(*S,sizeX1),
      X1b = MVT::Clone(*S,sizeX1),
      X2  = MVT::Clone(*S,sizeX2),
      Y2  = MVT::Clone(*S,sizeX2),
      X2b = MVT::Clone(*S,sizeX2);
    {
      int info, rank;
      MT err;
      SerialDenseMatrix<int,ST> yTx1(sizeX1,sizeX1), yTx2(sizeX2,sizeX2);
      Array<int> yTx1_piv(sizeX1), yTx2_piv(sizeX2);
      Teuchos::LAPACK<int,ST> lapack;
      TEUCHOS_TEST_FOR_EXCEPTION(sizeX1 < sizeX2,std::logic_error,"Internal logic error: sizeX1 < sizeX2.");
      Array<ST> LUwork(sizeX1);
      std::vector<MT> norms1(sizeX1), norms2(sizeX2);
      // use a BasicOrthoManager for testing
      RCP<MatOrthoManager<ST,MV,OP> > OM_basic = rcp( new BasicOrthoManager<ST,MV,OP>(M) );

      // Y1
      MVT::MvRandom(*Y1);
      rank = OM_basic->normalize(*Y1);
      TEUCHOS_TEST_FOR_EXCEPTION(rank != sizeX1, std::runtime_error,
          "normalize(Y1) returned rank " << rank << " from "
          << sizeX1 << " vectors. Cannot continue.");
      err = OM_basic->orthonormError(*Y1);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "normalize(Y1) did meet tolerance: orthonormError(Y1) == " << err);
      MyOM->stream(Warnings) << "   || <Y1,Y1> - I || : " << err << endl;

      // Y2
      MVT::MvRandom(*Y2);
      rank = OM_basic->normalize(*Y2);
      TEUCHOS_TEST_FOR_EXCEPTION(rank != sizeX2, std::runtime_error,
          "normalize(Y1) returned rank " << rank << " from "
          << sizeX2 << " vectors. Cannot continue.");
      err = OM_basic->orthonormError(*Y2);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "projectAndNormalize(Y2,Y1) did not meet tolerance: orthonormError(Y2) == " << err);
      MyOM->stream(Warnings) << "   || <Y2,Y2> - I || : " << err << endl;

      // X1 ortho to Y2
      MVT::MvRandom(*X1);
      OM_basic->project(*X1,tuple<RCP<const MV> >(Y2));
      err = OM_basic->orthogError(*X1,*Y2);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "project(X1,Y2) did not meet tolerance: orthog(X1,Y2) == " << err);
      MyOM->stream(Warnings) << "   || <X1,Y2> ||     : " << err << endl;
      MVT::MvNorm(*X1,norms1);
      for (unsigned int i=0; i<norms1.size(); i++) norms1[i] = ONE/norms1[i];
      MVT::MvScale(*X1,norms1);

      // Compute X1b so that <X1b,Y1> = I
      // Compute LU factorization of <Y1,X1> and use it to compute explicit inverse of <Y1,X1>
      OM_basic->innerProd(*Y1,*X1,yTx1);
      lapack.GETRF(yTx1.numRows(),yTx1.numCols(),yTx1.values(),yTx1.stride(),&yTx1_piv.front(),&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Error computing LU factorization of <Y1,X1>.");
      lapack.GETRI(yTx1.numRows(),yTx1.values(),yTx1.stride(),&yTx1_piv.front(),&LUwork.front(),LUwork.size(),&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Error computing LU inverse of <Y1,X1>.");
      // Set X1b=X1*inv(<Y1,X1>), so that <Y1,X1b> = <Y1,X1*inv(<Y1,X1>)> = <Y1,X1>*inv(<Y1,X1>) = I
      MVT::MvTimesMatAddMv(ONE,*X1,yTx1,ZERO,*X1b);
      // Test that this was successful
      OM_basic->innerProd(*Y1,*X1b,yTx1);
      for (int i=0; i<sizeX1; i++) {
        yTx1(i,i) -= ONE;
      }
      err = yTx1.normFrobenius();
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error, "Failed to make <X1b,Y1> = I.");
      // test <X1b,Y2> == 0
      err = OM_basic->orthogError(*X1b,*Y2);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "X1b did not meet tolerance: orthog(X1b,Y2) == " << err);
      MyOM->stream(Warnings) << "   || <X1b,Y2> ||     : " << err << endl;

      // Set X1=X1b*M1, so that <Y1,X1> is s.p.d.
      // Set M1 = R1'*R1, where R1 is random matrix
      SerialDenseMatrix<int,ST> R1(yTx1), M1(yTx1);
      Anasazi::randomSDM(R1);
      M1.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,R1,R1,ZERO);
      M1.scale(ONE/M1.normFrobenius());
      MVT::MvTimesMatAddMv(ONE,*X1b,M1,ZERO,*X1);
      // test <X1,Y2> == 0
      err = OM_basic->orthogError(*X1,*Y2);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "New X1 did not meet tolerance: orthog(X1,Y2) == " << err);
      MyOM->stream(Warnings) << "   || <X1,Y2> ||     : " << err << endl;
      // test explicit symmetry of <Y1,X1>:  yTx1 - yTx1' == 0
      SerialDenseMatrix<int,ST> xTy1(yTx1);
      OM_basic->innerProd(*Y1,*X1,yTx1);
      OM_basic->innerProd(*X1,*Y1,xTy1);
      xTy1 -= yTx1;
      err = xTy1.normFrobenius();
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "New <Y1,X1> not symmetric: ||<Y1,X1> - <X1,Y1>|| == " << err);
      MyOM->stream(Warnings) << "   ||<Y1,X1> - <X1,Y1>||     : " << err << endl;
      // test s.p.d. of <Y1,X1>: try to compute a cholesky
      lapack.POTRF('U',yTx1.numCols(),yTx1.values(),yTx1.stride(),&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::runtime_error,
          "New <Y1,X1> not s.p.d.: couldn't computed Cholesky: info == " << info
          << "\nyTx1: \n" << printMat(yTx1));

      // X2 ortho to Y1
      MVT::MvRandom(*X2);
      OM_basic->project(*X2,tuple<RCP<const MV> >(Y1));
      err = OM_basic->orthogError(*X2,*Y1);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "project(X2,Y1) did not meet tolerance: orthog(X2,Y1) == " << err);
      MyOM->stream(Warnings) << "   || <X2,Y1> ||     : " << err << endl;
      MVT::MvNorm(*X2,norms2);
      for (unsigned int i=0; i<norms2.size(); i++) norms2[i] = ONE/norms2[i];
      MVT::MvScale(*X2,norms2);

      // Compute X2b so that <X2b,Y2> = I
      // Compute LU factorization of <Y2,X2> and use it to compute explicit inverse of <Y2,X2>
      OM_basic->innerProd(*Y2,*X2,yTx2);
      lapack.GETRF(yTx2.numRows(),yTx2.numCols(),yTx2.values(),yTx2.stride(),&yTx2_piv.front(),&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Error computing LU factorization of <Y2,X2>.");
      lapack.GETRI(yTx2.numRows(),yTx2.values(),yTx2.stride(),&yTx2_piv.front(),&LUwork.front(),LUwork.size(),&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::logic_error, "Error computing LU inverse of <Y2,X2>.");
      // Set X2b=X2*inv(<Y2,X2>), so that <Y2,X2b> = <Y2,X2*inv(<Y2,X2>)> = <Y2,X2>*inv(<Y2,X2>) = I
      MVT::MvTimesMatAddMv(ONE,*X2,yTx2,ZERO,*X2b);
      // Test that this was successful
      OM_basic->innerProd(*Y2,*X2b,yTx2);
      for (int i=0; i<sizeX2; i++) {
        yTx2(i,i) -= ONE;
      }
      err = yTx2.normFrobenius();
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error, "Failed to make <X2b,Y2> = I.");
      // test <X2b,Y1> == 0
      err = OM_basic->orthogError(*X2b,*Y1);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "X2b did not meet tolerance: orthog(X2b,Y1) == " << err);
      MyOM->stream(Warnings) << "   || <X2b,Y1> ||     : " << err << endl;

      // Set X2=X2b*M2, so that <Y2,X2> is s.p.d.
      // Set M2 = R2'*R2, where R2 is random matrix
      SerialDenseMatrix<int,ST> R2(yTx2), M2(yTx2);
      Anasazi::randomSDM(R2);
      M2.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,R2,R2,ZERO);
      M2.scale(ONE/M2.normFrobenius());
      MVT::MvTimesMatAddMv(ONE,*X2b,M2,ZERO,*X2);
      // test <X2,Y1> == 0
      err = OM_basic->orthogError(*X2,*Y1);
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "New X2 did not meet tolerance: orthog(X2,Y1) == " << err);
      MyOM->stream(Warnings) << "   || <X2,Y1> ||     : " << err << endl;
      // test explicit symmetry of <Y2,X2>:  yTx2 - yTx2' == 0
      SerialDenseMatrix<int,ST> xTy2(yTx2);
      OM_basic->innerProd(*Y2,*X2,yTx2);
      OM_basic->innerProd(*X2,*Y2,xTy2);
      xTy2 -= yTx2;
      err = xTy2.normFrobenius();
      TEUCHOS_TEST_FOR_EXCEPTION(err > TOL,std::runtime_error,
          "New <Y2,X2> not symmetric: ||<Y2,X2> - <X2,Y2>|| == " << err);
      MyOM->stream(Warnings) << "   ||<Y2,X2> - <X2,Y2>||     : " << err << endl;
      // test s.p.d. of <Y2,X2>: try to compute a cholesky
      lapack.POTRF('U',yTx2.numCols(),yTx2.values(),yTx2.stride(),&info);
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,std::runtime_error,
          "New <Y2,X2> not s.p.d.: couldn't computed Cholesky: info == " << info
          << "\nyTx2: \n" << printMat(yTx2));
    }
    MyOM->stream(Warnings) << endl;


    {
      // just a random multivector
      MVT::MvRandom(*S);

      MyOM->stream(Errors) << " projectGen(): testing on random multivector " << endl;
      numFailed += testProjectGen(OM,S,X1,Y1,X2,Y2,false);
      MyOM->stream(Errors) << " projectGen(biOrtho): testing on random multivector " << endl;
      numFailed += testProjectGen(OM,S,X1b,Y1,X2b,Y2,true);
    }


    {
      // run a X1,Y2 range multivector against P_{X1,X1} P_{Y2,Y2}
      // note, this is allowed under the restrictions on projectGen,
      // because <X1,Y2> = 0
      // also, <Y2,Y2> = I, but <X1,X1> != I, so biOrtho must be set to false
      // it should destory the data, as:
      // P_{X1,X1} P_{Y2,Y2} (X1*C1 + Y2*C2) = P_{X1,X1} X1*C1 = 0
      // and
      // P_{Y2,Y2} P_{X1,X1} (X1*C1 + Y2*C2) = P_{Y2,Y2} Y2*C2 = 0
      SerialDenseMatrix<int,ST> C1(sizeX1,sizeS), C2(sizeX2,sizeS);
      Anasazi::randomSDM(C1);
      Anasazi::randomSDM(C2);
      MVT::MvTimesMatAddMv(ONE,*X1,C1,ZERO,*S);
      MVT::MvTimesMatAddMv(ONE,*Y2,C2,ONE,*S);

      MyOM->stream(Errors) << " projectGen(): testing [X1 Y2]-range multivector against P_{X1,X1} P_{Y2,Y2} " << endl;
      numFailed += testProjectGen(OM,S,X1,X1,Y2,Y2,false);
    }


    {
      // similar to above
      //
      // run a X2,Y1 range multivector against P_{X2,X2} P_{Y1,Y1}
      // note, this is allowed under the restrictions on projectGen,
      // because <X2,Y1> = 0
      // also, <Y1,Y1> = I, but <X2,X2> != I, so biOrtho must be set to false
      // it should destory the data, as:
      // P_{X2,X2} P_{Y1,Y1} (X2*C2 + Y1*C1) = P_{X2,X2} X2*C2 = 0
      // and
      // P_{Y1,Y1} P_{X2,X2} (X2*C2 + Y1*C1) = P_{Y1,Y1} Y1*C1 = 0
      SerialDenseMatrix<int,ST> C1(sizeX1,sizeS), C2(sizeX2,sizeS);
      Anasazi::randomSDM(C1);
      Anasazi::randomSDM(C2);
      MVT::MvTimesMatAddMv(ONE,*X2,C2,ZERO,*S);
      MVT::MvTimesMatAddMv(ONE,*Y1,C1,ONE,*S);

      MyOM->stream(Errors) << " projectGen(): testing [X2 Y1]-range multivector against P_{X2,X2} P_{Y1,Y1} " << endl;
      numFailed += testProjectGen(OM,S,X2,X2,Y1,Y1,false);
    }


    {
      // just a random multivector
      MVT::MvRandom(*S);

      MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on random multivector " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
      MyOM->stream(Errors) << " projectAndNormalizeGen(biOrtho): testing on random multivector " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
    }


    {
      // run a X1,Y2 range multivector against P_{X1,X1} P_{Y2,Y2}
      // note, this is allowed under the restrictions on projectAndNormalizeGen,
      // because <X1,Y2> = 0
      // also, <Y2,Y2> = I, but <X1,X1> != I, so biOrtho must be set to false
      // it should require randomization, as
      // P_{X1,X1} P_{Y2,Y2} (X1*C1 + Y2*C2) = P_{X1,X1} X1*C1 = 0
      // and
      // P_{Y2,Y2} P_{X1,X1} (X1*C1 + Y2*C2) = P_{Y2,Y2} Y2*C2 = 0
      SerialDenseMatrix<int,ST> C1(sizeX1,sizeS), C2(sizeX2,sizeS);
      Anasazi::randomSDM(C1);
      Anasazi::randomSDM(C2);
      MVT::MvTimesMatAddMv(ONE,*X1,C1,ZERO,*S);
      MVT::MvTimesMatAddMv(ONE,*Y2,C2,ONE,*S);

      MyOM->stream(Errors) << " projectAndNormalizeGen(): testing [X1 Y2]-range multivector against P_{X1,X1} P_{Y2,Y2} " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1,X1,Y2,Y2,false);
    }


    {
      // similar to above
      //
      // run a X2,Y1 range multivector against P_{X2,X2} P_{Y1,Y1}
      // note, this is allowed under the restrictions on projectGen,
      // because <X2,Y1> = 0
      // also, <Y1,Y1> = I, but <X2,X2> != I, so biOrtho must be set to false
      // it should require randomization, as
      // P_{X2,X2} P_{Y1,Y1} (X2*C2 + Y1*C1) = P_{X2,X2} X2*C2 = 0
      // and
      // P_{Y1,Y1} P_{X2,X2} (X2*C2 + Y1*C1) = P_{Y1,Y1} Y1*C1 = 0
      SerialDenseMatrix<int,ST> C1(sizeX1,sizeS), C2(sizeX2,sizeS);
      Anasazi::randomSDM(C1);
      Anasazi::randomSDM(C2);
      MVT::MvTimesMatAddMv(ONE,*X2,C2,ZERO,*S);
      MVT::MvTimesMatAddMv(ONE,*Y1,C1,ONE,*S);

      MyOM->stream(Errors) << " projectGen(): testing [X2 Y1]-range multivector against P_{X2,X2} P_{Y1,Y1} " << endl;
      numFailed += testProjectGen(OM,S,X2,X2,Y1,Y1,false);
    }


    if (sizeS > 2) {
      // rank-deficient
      MVT::MvRandom(*S);
      RCP<MV> mid = MVT::Clone(*S,1);
      SerialDenseMatrix<int,ST> c(sizeS,1);
      MVT::MvTimesMatAddMv(ONE,*S,c,ZERO,*mid);
      std::vector<int> ind(1);
      ind[0] = sizeS-1;
      MVT::SetBlock(*mid,ind,*S);

      MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on rank-deficient multivector " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
      MyOM->stream(Errors) << " projectAndNormalizeGen(biOrtho): testing on rank-deficient multivector " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
    }


    if (sizeS > 1) {
      // rank-1
      RCP<MV> one = MVT::Clone(*S,1);
      MVT::MvRandom(*one);
      SerialDenseMatrix<int,ST> scaleS(sizeS,1);
      Anasazi::randomSDM(scaleS);
      // put multiple of column 0 in columns 0:sizeS-1
      for (int i=0; i<sizeS; i++) {
        std::vector<int> ind(1);
        ind[0] = i;
        RCP<MV> Si = MVT::CloneViewNonConst(*S,ind);
        MVT::MvAddMv(scaleS(i,0),*one,ZERO,*one,*Si);
      }

      MyOM->stream(Errors) << " projectAndNormalizeGen(): testing on rank-1 multivector " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1,Y1,X2,Y2,false);
      MyOM->stream(Errors) << " projectAndNormalizeGen(biOrtho): testing on rank-1 multivector " << endl;
      numFailed += testProjectAndNormalizeGen(OM,S,X1b,Y1,X2b,Y2,true);
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,cout,success);

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

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
int testProjectAndNormalizeGen(RCP<GenOrthoManager<ST,MV,OP> > OM,
                               RCP<const MV> S,
                               RCP<const MV> X1, RCP<const MV> Y1,
                               RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho)
{
  typedef MultiVecTraits<double,MV>    MVT;
  typedef OperatorTraits<double,MV,OP> OPT;
  typedef ScalarTraits<double>         SCT;
  typedef SCT::magnitudeType           MT;

  const ST ONE = SCT::one();
  const MT ZERO = SCT::magnitude(SCT::zero());
  const int sizeS = MVT::GetNumberVecs(*S);
  const int sizeX1 = MVT::GetNumberVecs(*X1);
  const int sizeX2 = MVT::GetNumberVecs(*X2);
  int numerr = 0;
  bool hasM = (OM->getOp() != null);
  std::ostringstream sout;

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
  // if hasM:
  // with and without MX1,MY1   (1---)
  // with and without MX2,MY2  (1----)
  // with and without MS      (1-----)
  //
  // as hasM controls the upper level bits, we need only run test cases 0-7 if hasM==false
  // otherwise, we run test cases 0-63
  //

  int numtests;
  RCP<const MV> MX1, MX2, MY1, MY2;
  RCP<MV> MS;
  if (hasM) {
    MX1 = MVT::Clone(*S,sizeX1);
    MY1 = MVT::Clone(*S,sizeX1);
    MX2 = MVT::Clone(*S,sizeX2);
    MY2 = MVT::Clone(*S,sizeX2);
    MS  = MVT::Clone(*S,sizeS);
    OPT::Apply(*(OM->getOp()),*X1,const_cast<MV&>(*MX1));
    OPT::Apply(*(OM->getOp()),*Y1,const_cast<MV&>(*MY1));
    OPT::Apply(*(OM->getOp()),*X2,const_cast<MV&>(*MX2));
    OPT::Apply(*(OM->getOp()),*Y2,const_cast<MV&>(*MY2));
    OPT::Apply(*(OM->getOp()),*S ,*MS);
    numtests = 64;
  }
  else {
    numtests = 8;
  }

  // test ortho error before orthonormalizing
  if (Y1 != null) {
    MT err = OM->orthogError(*S,*Y1);
    sout << "   || <S,Y1> || before     : " << err << endl;
  }
  if (Y2 != null) {
    MT err = OM->orthogError(*S,*Y2);
    sout << "   || <S,Y2> || before     : " << err << endl;
  }

  for (int t=0; t<numtests; t++) {

    // pointers to simplify calls below
    RCP<const MV> lclMX1;
    RCP<const MV> lclMY1;
    RCP<const MV> lclMX2;
    RCP<const MV> lclMY2;
    RCP<MV> lclMS;
    if ( t & 8 ) {
      lclMX1 = MX1;
      lclMY1 = MY1;
    }
    if ( t & 16 ) {
      lclMX2 = MX2;
      lclMY2 = MY2;
    }
    if ( t & 32 ) {
      lclMS = MS;
    }

    Array<RCP<const MV> > theX, theY, theMX, theMY;
    RCP<SerialDenseMatrix<int,ST> > B = rcp( new SerialDenseMatrix<int,ST>(sizeS,sizeS) );
    Array<RCP<SerialDenseMatrix<int,ST> > > C;
    if ( (t && 3) == 0 ) {
      // neither <X1,Y1> nor <X2,Y2>
      // C, theX and theY are already empty
    }
    else if ( (t && 3) == 1 ) {
      // <X1,Y1>
      theX = tuple(X1);
      theY = tuple(Y1);
      theMX = tuple(MX1);
      theMY = tuple(MY1);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)) );
    }
    else {
      // <X1,Y2> and <X2,Y2>, and the reverse.
      theX = tuple(X1,X2);
      theY = tuple(Y1,Y2);
      theMX = tuple(MX1,MX2);
      theMY = tuple(MY1,MY2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)),
                 rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }

    bool localIsBiortho = isBiortho;
    if (localIsBiortho) {
      if ( t && 4 == 0 ) {
        localIsBiortho = false;
      }
    }

    try {
      // call routine
      // if (t && 3) == 3, {
      //    call with reversed input: <X2,Y2> <X1,Y1>
      // }
      // test all outputs for correctness
      // test all outputs for equivalence

      // here is where the outputs go
      Array<RCP<const MV> > S_outs;
      Array<RCP<const MV> > MS_outs;
      Array<Array<RCP<SerialDenseMatrix<int,ST> > > > C_outs;
      Array<RCP<SerialDenseMatrix<int,ST> > > B_outs;
      RCP<MV> Scopy, MScopy;
      Array<int> ret_out;

      // copies of S,MS
      Scopy = MVT::CloneCopy(*S);
      if (lclMS != Teuchos::null) {
        MScopy = MVT::CloneCopy(*lclMS);
      }
      // randomize this data, it should be overwritten
      Anasazi::randomSDM(*B);
      for (unsigned int i=0; i<C.size(); i++) {
        Anasazi::randomSDM(*C[i]);
      }
      // run test
      int ret = OM->projectAndNormalizeGen(
            *Scopy,theX,theY,localIsBiortho,C,B,MScopy,theMX,theMY);
      sout << "projectAndNormalizeGen() returned rank " << ret << endl;
      if (ret == 0) {
        sout << "   Cannot continue." << endl;
        numerr++;
        break;
      }
      ret_out.push_back(ret);
      // projectAndNormalizeGen() is only required to return a
      // basis of rank "ret"
      // this is what we will test:
      //   the first "ret" columns in Scopy, MScopy
      //   the first "ret" rows in B
      // save just the parts that we want
      // we allocate S and MS for each test, so we can save these as views
      // however, save copies of the C and B
      if (ret < sizeS) {
        std::vector<int> ind(ret);
        for (int i=0; i<ret; i++) {
          ind[i] = i;
        }
        S_outs.push_back( MVT::CloneView(*Scopy,ind) );
        if (MScopy != null) {
          MS_outs.push_back( MVT::CloneView(*MScopy,ind) );
        }
        else {
          MS_outs.push_back( Teuchos::null );
        }
        B_outs.push_back( rcp( new SerialDenseMatrix<int,ST>(Teuchos::Copy,*B,ret,sizeS) ) );
      }
      else {
        S_outs.push_back( Scopy );
        MS_outs.push_back( MScopy );
        B_outs.push_back( rcp( new SerialDenseMatrix<int,ST>(*B) ) );
      }
      C_outs.push_back( Array<RCP<SerialDenseMatrix<int,ST> > >(0) );
      if (C.size() > 0) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[0]) ) );
      }
      if (C.size() > 1) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[1]) ) );
      }

      // test all outputs for correctness
      for (unsigned int o=0; o<S_outs.size(); o++) {
        // MS == M*S
        if (MS_outs[o] != null) {
          RCP<MV> tmp = MVT::Clone(*S_outs[o],ret_out[o]);
          OPT::Apply(*(OM->getOp()),*S_outs[o],*tmp);
          MT err = MVDiff(*tmp,*MS_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| MS - M*S ||           : " << err << endl;
        }
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
        // <Y1,S> == 0
        if (theY.size() > 0 && theY[0] != null) {
          MT err = OM->orthogError(*theY[0],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <Y[0],S> || after      : " << err << endl;
        }
        // <Y2,S> == 0
        if (theY.size() > 1 && theY[1] != null) {
          MT err = OM->orthogError(*theY[1],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <Y[1],S> || after      : " << err << endl;
        }
      }

      // test all outputs for consistency
      // finish

    }
    catch (const OrthoError &e) {
      sout << "   -------------------------------------------         projectAndNormalizeGen() threw exception" << endl;
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
int testProjectGen(RCP<GenOrthoManager<ST,MV,OP> > OM,
                   RCP<const MV> S,
                   RCP<const MV> X1, RCP<const MV> Y1,
                   RCP<const MV> X2, RCP<const MV> Y2, bool isBiortho)
{
  typedef MultiVecTraits<double,MV>    MVT;
  typedef OperatorTraits<double,MV,OP> OPT;
  typedef ScalarTraits<double>         SCT;
  typedef SCT::magnitudeType           MT;

  const ST ONE = SCT::one();
  const int sizeS = MVT::GetNumberVecs(*S);
  const int sizeX1 = MVT::GetNumberVecs(*X1);
  const int sizeX2 = MVT::GetNumberVecs(*X2);
  int numerr = 0;
  bool hasM = (OM->getOp() != null);
  std::ostringstream sout;

  //
  // output tests:
  //   <S_out,Y1> = 0
  //   <S_out,Y2> = 0
  //   S_in = S_out + X1 C1 + X2 C2
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
  // with C1 and C2
  // with and without isBiortho (if isBiortho==true)  (1--)
  //
  // if hasM:
  // with and without MX1,MY1   (1---)
  // with and without MX2,MY2  (1----)
  // with and without MS      (1-----)
  //
  // as hasM controls the upper level bits, we need only run test cases 0-7 if hasM==false
  // otherwise, we run test cases 0-63
  //

  int numtests;
  RCP<const MV> MX1, MX2, MY1, MY2;
  RCP<MV> MS;
  if (hasM) {
    MX1 = MVT::Clone(*S,sizeX1);
    MY1 = MVT::Clone(*S,sizeX1);
    MX2 = MVT::Clone(*S,sizeX2);
    MY2 = MVT::Clone(*S,sizeX2);
    MS  = MVT::Clone(*S,sizeS);
    OPT::Apply(*(OM->getOp()),*X1,const_cast<MV&>(*MX1));
    OPT::Apply(*(OM->getOp()),*Y1,const_cast<MV&>(*MY1));
    OPT::Apply(*(OM->getOp()),*X2,const_cast<MV&>(*MX2));
    OPT::Apply(*(OM->getOp()),*Y2,const_cast<MV&>(*MY2));
    OPT::Apply(*(OM->getOp()),*S ,*MS);
    numtests = 64;
  }
  else {
    numtests = 8;
  }

  // test ortho error before orthonormalizing
  if (Y1 != null) {
    MT err = OM->orthogError(*S,*Y1);
    sout << "   || <S,Y1> || before     : " << err << endl;
  }
  if (Y2 != null) {
    MT err = OM->orthogError(*S,*Y2);
    sout << "   || <S,Y2> || before     : " << err << endl;
  }

  for (int t=0; t<numtests; t++) {

    // pointers to simplify calls below
    RCP<const MV> lclMX1;
    RCP<const MV> lclMY1;
    RCP<const MV> lclMX2;
    RCP<const MV> lclMY2;
    RCP<MV> lclMS;
    if ( t & 8 ) {
      lclMX1 = MX1;
      lclMY1 = MY1;
    }
    if ( t & 16 ) {
      lclMX2 = MX2;
      lclMY2 = MY2;
    }
    if ( t & 32 ) {
      lclMS = MS;
    }

    Array<RCP<const MV> > theX, theY, theMX, theMY;
    Array<RCP<SerialDenseMatrix<int,ST> > > C;
    if ( (t && 3) == 0 ) {
      // neither <X1,Y1> nor <X2,Y2>
      // C, theX and theY are already empty
    }
    else if ( (t && 3) == 1 ) {
      // <X1,Y1>
      theX = tuple(X1);
      theY = tuple(Y1);
      theMX = tuple(MX1);
      theMY = tuple(MY1);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)) );
    }
    else {
      // <X1,Y2> and <X2,Y2>, and the reverse.
      theX = tuple(X1,X2);
      theY = tuple(Y1,Y2);
      theMX = tuple(MX1,MX2);
      theMY = tuple(MY1,MY2);
      C = tuple( rcp(new SerialDenseMatrix<int,ST>(sizeX1,sizeS)),
                 rcp(new SerialDenseMatrix<int,ST>(sizeX2,sizeS)) );
    }

    bool localIsBiortho = isBiortho;
    if (localIsBiortho) {
      if ( t && 4 == 0 ) {
        localIsBiortho = false;
      }
    }

    try {
      // call routine
      // if (t && 3) == 3, {
      //    call with reversed input: <X2,Y2> <X1,Y1>
      // }
      // test all outputs for correctness
      // test all outputs for equivalence

      // here is where the outputs go
      Array<RCP<MV> > S_outs;
      Array<RCP<MV> > MS_outs;
      Array<Array<RCP<SerialDenseMatrix<int,ST> > > > C_outs;
      RCP<MV> Scopy, MScopy;

      // copies of S,MS
      Scopy = MVT::CloneCopy(*S);
      if (lclMS != Teuchos::null) {
        MScopy = MVT::CloneCopy(*lclMS);
      }
      // randomize this data, it should be overwritten
      for (unsigned int i=0; i<C.size(); i++) {
        Anasazi::randomSDM(*C[i]);
      }
      // run test
      OM->projectGen(*Scopy,theX,theY,localIsBiortho,C,MScopy,theMX,theMY);
      // we allocate S and MS for each test, so we can save these as views
      // however, save copies of the C
      S_outs.push_back( Scopy );
      MS_outs.push_back( MScopy );
      C_outs.push_back( Array<RCP<SerialDenseMatrix<int,ST> > >(0) );
      if (C.size() > 0) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[0]) ) );
      }
      if (C.size() > 1) {
        C_outs.back().push_back( rcp( new SerialDenseMatrix<int,ST>(*C[1]) ) );
      }

      // test all outputs for correctness
      for (unsigned int o=0; o<S_outs.size(); o++) {
        // MS == M*S
        if (MS_outs[o] != null) {
          RCP<MV> tmp = MVT::Clone(*S_outs[o],sizeS);
          OPT::Apply(*(OM->getOp()),*S_outs[o],*tmp);
          MT err = MVDiff(*tmp,*MS_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| MS - M*S ||           : " << err << endl;
        }
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
        // <Y1,S> == 0
        if (theY.size() > 0 && theY[0] != null) {
          MT err = OM->orthogError(*theY[0],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <Y[0],S> || after      : " << err << endl;
        }
        // <Y2,S> == 0
        if (theY.size() > 1 && theY[1] != null) {
          MT err = OM->orthogError(*theY[1],*S_outs[o]);
          if (err > TOL) {
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
          sout << "  " << t << "|| <Y[1],S> || after      : " << err << endl;
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
            sout << "         vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
            numerr++;
          }
        }
      }

    }
    catch (const OrthoError &e) {
      sout << "   -------------------------------------------         projectGen() threw exception" << endl;
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



double MVDiff(const MV &X, const MV &Y)
{
  typedef MultiVecTraits<double,MV>    MVT;
  //typedef OperatorTraits<double,MV,OP> OPT; // unused
  typedef ScalarTraits<double>         SCT;
  typedef SCT::magnitudeType           MT;

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
