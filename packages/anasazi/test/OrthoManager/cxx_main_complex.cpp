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
//  This test is for the BasicOrthoManager
//
// The matrix used is from MatrixMarket:
// Name: MHD1280B: Alfven Spectra in Magnetohydrodynamics
// Source: Source: A. Booten, M.N. Kooper, H.A. van der Vorst, S. Poedts and J.P. Goedbloed University of Utrecht, the Netherlands
// Discipline: Plasma physics
// URL: http://math.nist.gov/MatrixMarket/data/NEP/mhd/mhd1280b.html
// Size: 1280 x 1280
// NNZ: 22778 entries

#include "AnasaziConfigDefs.hpp"
#include "AnasaziModalSolverUtils.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "iohb.h"
#endif

// templated multivector and sparse matrix classes
#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"


#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace Teuchos;
using namespace Anasazi;

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;
#else
  typedef double ST;
  // no complex. quit with failure.
int main() {
  cout << "Not compiled with complex support." << endl;
  cout << "End Result: TEST FAILED" << endl;
  return -1;
}
#endif
typedef ScalarTraits<ST>         SCT;
typedef SCT::magnitudeType        MT;
typedef MultiVec<ST>              MV;
typedef Operator<ST>              OP;
typedef MultiVecTraits<ST,MV>    MVT;
typedef OperatorTraits<ST,MV,OP> OPT;

const ST ONE = SCT::one();
const MT ZERO = SCT::magnitude(SCT::zero());

// this is the tolerance that all tests are performed against
const MT TOL = 1.0e-13;
const MT ATOL = 100;

// declare an output manager for handling local output
RefCountPtr< Anasazi::BasicOutputManager<ST> > MyOM;

// some forward declarations
int testProjectAndNormalize(RefCountPtr<MatOrthoManager<ST,MV,OP> > OM, 
                            RefCountPtr<MV> X, RefCountPtr<const MV> Q);
int testProject(RefCountPtr<MatOrthoManager<ST,MV,OP> > OM, 
                RefCountPtr<const MV> X, RefCountPtr<const MV> Q);
int testNormalize(RefCountPtr<MatOrthoManager<ST,MV,OP> > OM, RefCountPtr<MV> X);
MT MVDiff(const MV &X, const MV &Y);

int main(int argc, char *argv[]) 
{
  int info = 0;
  int MyPID;
  
#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &MyPID );
#else
  MyPID = 0;
#endif

  bool verbose = false;
  int numFailed = 0;
  bool debug = false;
  std::string filename("mhd1280b.cua");
  std::string which("LR");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SR or LR).");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

#ifndef HAVE_ANASAZI_TRIUTILS
  cout << "This test requires Triutils. Please configure with --enable-triutils." << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  cout << "End Result: TEST FAILED" << endl;	
  return -1;
#else

  // instantiate the output manager
  MyOM = rcp( new BasicOutputManager<ST>() );
  if (verbose) {
    // output in this driver will be sent to Anasazi::Warnings
    MyOM->setVerbosity(Anasazi::Warnings);
  }

  // Output Anasazi version
  MyOM->stream(Anasazi::Warnings) << Anasazi_Version() << endl << endl;

  // Get the data from the HB file
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;
  info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,
                              &colptr,&rowind,&dvals);
  if (info == 0 || nnz < 0) {
    cout << "Error reading '" << filename << "'" << endl;
    cout << "End Result: TEST FAILED" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  // Convert interleaved doubles to complex values
  std::vector<ST> cvals(nnz);
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ST(dvals[ii*2],dvals[ii*2+1]);
  }
  // Build the problem matrix
  RefCountPtr< MyBetterOperator<ST> > M 
    = rcp( new MyBetterOperator<ST>(dim,colptr,nnz,rowind,&cvals[0]) );

  const int sizeX = 4;
  const int sizeQ = 10; // MUST: sizeX + sizeQ <= elements[0]-1
  MT err;
  bool tfail;


  // Create an SVQB ortho manager 
  RefCountPtr<MatOrthoManager<ST,MV,OP> > OM   = rcp( new SVQBOrthoManager<ST,MV,OP>(Teuchos::null,debug) );
  RefCountPtr<MatOrthoManager<ST,MV,OP> > OM_M = rcp( new SVQBOrthoManager<ST,MV,OP>(M,debug) );

  // multivector to spawn off of
  RefCountPtr<MV> X = rcp( new MyMultiVec<ST>(dim, sizeX) );

  MyOM->stream(Errors) << " Generating Q1 for project() : testing... " << endl;
  RefCountPtr<MV> Q1  = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q1);
  int dummy = 0;
  dummy = OM_M->normalize(*Q1,null);
  tfail = false;
  if ( dummy != sizeQ ) {
    MyOM->stream(Errors) << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! normalize() returned " << dummy << " . Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  err = OM_M->orthonormError(*Q1);
  if (err > TOL) {
    MyOM->stream(Errors) << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  if (tfail) {
    MyOM->stream(Warnings) << "   || Q1^T M Q1 - I ||_F        : " << err << endl << endl;
  }

  MyOM->stream(Errors) << " Generating Q2 for project() : testing... " << endl;
  RefCountPtr<MV> Q2 = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q2);
  tfail = false;
  dummy = OM->normalize(*Q2,null);
  if ( dummy != sizeQ ) {
    MyOM->stream(Errors) << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! normalize() returned " << dummy << " . Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  err = OM->orthonormError(*Q2);
  if (err > TOL) {
    MyOM->stream(Errors) << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  if (tfail) {
    MyOM->stream(Warnings) << "   || Q2^T Q2 - I ||_F          : " << err << endl << endl;
  }

  {
    MyOM->stream(Errors) << " project() with M    : testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(OM_M,X,Q1);
  }
  {
    MyOM->stream(Errors) << " project() without M : testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(OM,X,Q2);
  }

  {
    MyOM->stream(Errors) << " normalize() with M    : testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(OM_M,X);
  }
  {
    MyOM->stream(Errors) << " normalize() without M : testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(OM,X);
  }


  {
    MyOM->stream(Errors) << " normalize() with M    : testing on rank-deficient multivector " << endl;
    std::vector<int> ind(1); 
    // Assuming here that sizeX == 4
    MVT::MvRandom(*X);
    ind[0] = 1;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    ind[0] = 2;
    MVT::SetBlock(*mid,ind,*X);
    
    numFailed += testNormalize(OM_M,X);
  }
  {
    MyOM->stream(Errors) << " normalize() without M : testing on rank-deficient multivector " << endl;
    std::vector<int> ind(1); 
    // Assuming here that sizeX == 4
    MVT::MvRandom(*X);
    ind[0] = 1;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    ind[0] = 2;
    MVT::SetBlock(*mid,ind,*X);
    
    numFailed += testNormalize(OM,X);
  }


  {
    MyOM->stream(Errors) << " normalize() with M    : testing on rank-1 multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    // get column 0
    ind[0] = 0;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    // put column 0 in columns 1:sizeX-1
    for (int i=1; i<sizeX; i++) {
      ind[0] = i;
      MVT::SetBlock(*mid,ind,*X);
    }
    
    numFailed += testNormalize(OM_M,X);
  }
  {
    MyOM->stream(Errors) << " normalize() without M : testing on rank-1 multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    // get column 0
    ind[0] = 0;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    // put column 0 in columns 1:sizeX-1
    for (int i=1; i<sizeX; i++) {
      ind[0] = i;
      MVT::SetBlock(*mid,ind,*X);
    }
    
    numFailed += testNormalize(OM,X);
  }

  {
    MyOM->stream(Errors) << " projectAndNormalize() with M    : testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(OM_M,X,Q1);
  }
  {
    MyOM->stream(Errors) << " projectAndNormalize() without M : testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(OM,X,Q2);
  }

  {
    MyOM->stream(Errors) << " projectAndNormalize() with M    : testing on Q-range multivector " << endl;
    SerialDenseMatrix<int,ST> B(sizeQ,sizeX);
    B.random();
    MVT::MvTimesMatAddMv(ONE,*Q1,B,ZERO,*X);
    
    numFailed += testProjectAndNormalize(OM_M,X,Q1);
  }
  {
    MyOM->stream(Errors) << " projectAndNormalize() without M : testing on Q-range multivector " << endl;
    SerialDenseMatrix<int,ST> B(sizeQ,sizeX);
    B.random();
    MVT::MvTimesMatAddMv(ONE,*Q2,B,ZERO,*X);
    
    numFailed += testProjectAndNormalize(OM,X,Q2);
  }

  {
    MyOM->stream(Errors) << " projectAndNormalize() with M    : testing on rank-deficient multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    // Assuming here that sizeX == 4
    ind[0] = 1;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    ind[0] = 2;
    MVT::SetBlock(*mid,ind,*X);
    
    numFailed += testProjectAndNormalize(OM_M,X,Q1);
  }
  {
    MyOM->stream(Errors) << " projectAndNormalize() without M : testing on rank-deficient multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    // Assuming here that sizeX == 4
    ind[0] = 1;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    ind[0] = 2;
    MVT::SetBlock(*mid,ind,*X);
    
    numFailed += testProjectAndNormalize(OM,X,Q2);
  }

  {
    MyOM->stream(Errors) << " projectAndNormalize() with M    : testing on rank-1 multivector " << endl;
    MVT::MvRandom(*X);
    std::vector<int> ind(1); 
    // get column 0
    ind[0] = 0;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    // put column 0 in columns 1:sizeX-1
    for (int i=1; i<sizeX; i++) {
      ind[0] = i;
      MVT::SetBlock(*mid,ind,*X);
    }
    
    numFailed += testProjectAndNormalize(OM_M,X,Q1);
  }
  {
    MyOM->stream(Errors) << " projectAndNormalize() without M : testing on rank-1 multivector " << endl;
    MVT::MvRandom(*X);
    std::vector<int> ind(1); 
    // get column 0
    ind[0] = 0;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    // put column 0 in columns 1:sizeX-1
    for (int i=1; i<sizeX; i++) {
      ind[0] = i;
      MVT::SetBlock(*mid,ind,*X);
    }
    
    numFailed += testProjectAndNormalize(OM,X,Q2);
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  // Clean up.
  free( dvals );
  free( colptr );
  free( rowind );

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
#endif  // isdef ANASAZI_TRIUTILS
}	




//////////////////////////////////////////////////////////////////////
int testProjectAndNormalize(RefCountPtr<MatOrthoManager<ST,MV,OP> > OM, 
                            RefCountPtr<MV> X, RefCountPtr<const MV> Q) {

  const int sizeX = MVT::GetNumberVecs(*X);
  const int sizeQ = MVT::GetNumberVecs(*Q);
  MT err;
  int rank;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  int numerr = 0;
  bool hasM = (OM->getOp() != null);
  std::ostringstream sout;
  RefCountPtr<MV> xcopy, mxcopy, smlOX;
  RefCountPtr<SerialDenseMatrix<int,ST> > C, smlC, R, smlR, newR;
  bool warning = false;

  int numtests;
  RefCountPtr<MV> MX;
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
  err = OM->orthonormError(*X);
  sout << "   || X^T M X - I ||_F before : " << err << endl;
  err = OM->orthogError(*X,*Q);
  sout << "   || Q^T M X ||_F before     : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<ST,MV,OP> MSU(rcp(new BasicOutputManager<ST>()));
  xcopy  = MVT::CloneCopy(*X);
  mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,OM->getOp().get(),*Q,sizeX,0);
  sout << "   MSU.massOrthonormalize returned " << iret << endl;
  sout << "   MSU.massOrthonormalize || X^T M X - I ||_F : " << OM->orthonormError(*xcopy) << endl;
  sout << "   MSU.massOrthonormalize || Q^T M X ||_F     : " << OM->orthogError(*Q,*xcopy) << endl;

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
      C->putScalar(ZERO);
      R = rcp( new SerialDenseMatrix<int,ST>(sizeX,sizeX) );
      R->putScalar(ZERO);
    }
    else {
      C = null;
      R = null;
    }

    try {
      rank = OM->projectAndNormalize(*xcopy,mxcopy,
                                     tuple< RefCountPtr< SerialDenseMatrix<int,ST> > >(C),R,
                                     tuple<RefCountPtr<const MV> >(Q) );

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
          warning = true;
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! warning!" << endl;
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
          warning = true;
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! warning!" << endl;
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
  if (numerr>0 || warning) type = Errors;
  MyOM->stream(type) << sout.str();
  MyOM->stream(type) << endl;

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testProject(RefCountPtr<MatOrthoManager<ST,MV,OP> > OM, 
                RefCountPtr<const MV> X, RefCountPtr<const MV> Q) {

  const int sizeX = MVT::GetNumberVecs(*X);
  const int sizeQ = MVT::GetNumberVecs(*Q);
  MT err;
  RefCountPtr<MV> tmp;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;
  std::ostringstream sout;
  RefCountPtr<MV> xcopy, mxcopy;
  RefCountPtr<SerialDenseMatrix<int,ST> > C;
  bool warning = false;

  int numtests;
  RefCountPtr<MV> MX;
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

  // first, run with MSUtils
  ModalSolverUtils<ST,MV,OP> MSU(rcp(new BasicOutputManager<ST>()));
  xcopy = MVT::CloneCopy(*X);
  mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,OM->getOp().get(),*Q,sizeX,1);
  sout << "   MSU.massOrthonormalize returned " << iret << endl;
  sout << "   MSU.massOrthonormalize error: " << OM->orthogError(*Q,*xcopy) << endl;

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
    }
    else {
      C = null;
    }

    try {
      OM->project(*xcopy,mxcopy,
                  tuple<RefCountPtr<SerialDenseMatrix<int,ST> > >(C),
                  tuple<RefCountPtr<const MV> >(Q));
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



//////////////////////////////////////////////////////////////////////
int testNormalize(RefCountPtr<MatOrthoManager<ST,MV,OP> > OM, RefCountPtr<MV> X) {

  const int sizeX = MVT::GetNumberVecs(*X);
  MT err;
  int rank;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;
  std::ostringstream sout;
  RefCountPtr<MV> xcopy, mxcopy, smlOX;
  RefCountPtr<SerialDenseMatrix<int,ST> > R, smlR, newR;
  bool warning = false;

  int numtests;
  RefCountPtr<MV> MX;
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
  err = OM->orthonormError(*X);
  sout << "   || X^T M X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<ST,MV,OP> MSU(rcp(new BasicOutputManager<ST>()));
  xcopy  = MVT::CloneCopy(*X);
  mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,OM->getOp().get(),*xcopy,sizeX,2);
  sout << "   MSU.massOrthonormalize returned " << iret << endl;
  sout << "   MSU.massOrthonormalize error: " << OM->orthonormError(*xcopy) << endl;

  for (int t=0; t<numtests; t++) {

    xcopy = MVT::CloneCopy(*X);

    if (t & 2) {
      mxcopy = MVT::CloneCopy(*MX);
    }
    else {
      mxcopy = null;
    }

    if (t & 1) {
      R = rcp( new SerialDenseMatrix<int,ST>(sizeX,sizeX) );
      R->putScalar(ZERO);
    }
    else {
      R = null;
    }

    try {
      rank = OM->normalize(*xcopy,mxcopy,R);
      sout << "normalize() returned rank " << rank << endl;
  
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
          warning = true;
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! warning!" << endl;
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
    }
    catch (OrthoError e) {
      sout << "   -------------------------------------------         normalize() Failed" << endl;
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
  RefCountPtr<MV> tmp = MVT::CloneCopy(X);

  // tmp <- tmp - Y
  MVT::MvAddMv(-ONE,Y,ONE,*tmp,*tmp);

  MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
  MT err = 0;
  for (int i=0; i<sizeX; i++) {
    err += SCT::magnitude(xTmx(i,i));
  }
  return SCT::magnitude(SCT::squareroot(err));
}
