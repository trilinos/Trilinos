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
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziModalSolverUtils.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

using namespace Teuchos;
using namespace Anasazi;

typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;
typedef MultiVecTraits<double,MV>    MVT;
typedef OperatorTraits<double,MV,OP> OPT;
typedef ScalarTraits<double>         SCT;
typedef SCT::magnitudeType   MT;

const double ONE = SCT::one();
const double ZERO = SCT::magnitude(SCT::zero());

// this is the tolerance that all tests are performed against
const double TOL = 1.0e-14;

bool verbose = false;


// some forward declarations
int testProjectAndNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                            RefCountPtr<MV> X, RefCountPtr<const MV> Q);
int testProject(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                RefCountPtr<const MV> X, RefCountPtr<const MV> Q);
int testNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X);
double MVDiff(const MV &X, const MV &Y);

int main(int argc, char *argv[]) 
{
  int i;
  int info = 0;
  
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  int MyPID = Comm.MyPID();
  
  int numFailed = 0;
  verbose = false;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }
  verbose = verbose && (MyPID==0);

  if (verbose) {
    cout << Anasazi_Version() << endl << endl;
  }
  
  ReturnType ret;
  
  // Problem information
  const int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;
  const int sizeX = 4;
  const int sizeQ = 10; // MUST: sizeX + sizeQ <= elements[0]-1
  double err;
  bool tfail;

  // Create problem
  RefCountPtr<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  // Get the mass matrix
  RefCountPtr<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  // Create a basic ortho manager 
  RefCountPtr<MatOrthoManager<double,MV,OP> > OM   = rcp( new BasicOrthoManager<double,MV,OP>() );
  RefCountPtr<MatOrthoManager<double,MV,OP> > OM_M = rcp( new BasicOrthoManager<double,MV,OP>(M) );

  // multivector to spawn off of
  RefCountPtr<MV> X = rcp( new Epetra_MultiVector(M->OperatorDomainMap(), sizeX) );

  cout << " Generating Q1 for project() testing... " << endl;
  RefCountPtr<MV> Q1  = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q1);
  int dummy = 0;
  ret = OM_M->normalize(*Q1,null,dummy);
  tfail = false;
  if ( ret  != Ok ) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! normalize() return " << ret << "/" << dummy << " . Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  err = OM_M->orthonormError(*Q1);
  if (err > TOL) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  if (verbose || tfail) {
    cout << "   || Q1^T M Q1 - I ||_F        : " << err << endl << endl;
  }

  cout << " Generating Q2 for project() testing... " << endl;
  RefCountPtr<MV> Q2 = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q2);
  tfail = false;
  ret = OM->normalize(*Q2,null,dummy);
  if ( ret != Ok ) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! normalize() return " << ret << "/" << dummy << " . Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  err = OM->orthonormError(*Q2);
  if (err > TOL) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
    tfail = true;
  }
  if (verbose || tfail) {
    cout << "   || Q2^T Q2 - I ||_F          : " << err << endl << endl;
  }

  {
    cout << " project() with M testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(OM_M,X,Q1);
  }
  {
    cout << " project() without M testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(OM,X,Q2);
  }

  {
    cout << " normalize() with M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(OM_M,X);
  }
  {
    cout << " normalize() without M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(OM,X);
  }


  {
    cout << " normalize() with M testing on rank-deficient multivector " << endl;
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
    cout << " normalize() without M testing on rank-deficient multivector " << endl;
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
    cout << " normalize() with M testing on rank-1 multivector " << endl;
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
    cout << " normalize() without M testing on rank-1 multivector " << endl;
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
    cout << " projectAndNormalize() with M testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(OM_M,X,Q1);
  }
  {
    cout << " projectAndNormalize() without M testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(OM,X,Q2);
  }

  {
    cout << " projectAndNormalize() with M testing on rank-deficient multivector " << endl;
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
    cout << " projectAndNormalize() without M testing on rank-deficient multivector " << endl;
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
    cout << " projectAndNormalize() with M testing on rank-1 multivector " << endl;
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
    cout << " projectAndNormalize() without M testing on rank-1 multivector " << endl;
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

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (numFailed) {
    cout << numFailed << " errors." << endl;
    cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  cout << "End Result: TEST PASSED" << endl;
  return 0;
}	




//////////////////////////////////////////////////////////////////////
int testProjectAndNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                            RefCountPtr<MV> X, RefCountPtr<const MV> Q) {

  const int sizeX = MVT::GetNumberVecs(*X);
  const int sizeQ = MVT::GetNumberVecs(*Q);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  int numerr = 0;
  bool hasM = (OM->getOp() != null);
  std::ostringstream sout;
  RefCountPtr<MV> xcopy, mxcopy, smlOX;
  RefCountPtr<SerialDenseMatrix<int,double> > C, smlC, R, smlR, newR;

  int numtests;
  RefCountPtr<MV> MX;
  if (hasM) {
    numtests = 4;
    MX = MVT::Clone(*X,sizeX);
    if ( OPT::Apply(*(OM->getOp()),*X,*MX) != Ok ) return 1;
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
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
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
      C = rcp( new SerialDenseMatrix<int,double>(sizeQ,sizeX) );
      C->putScalar(ZERO);
      R = rcp( new SerialDenseMatrix<int,double>(sizeX,sizeX) );
      R->putScalar(ZERO);
    }
    else {
      C = null;
      R = null;
    }

    ret = OM->projectAndNormalize(*xcopy,mxcopy,C,R,*Q,rank);
    switch (ret) {
    case Ok:
      // check the following:
      // rank == sizeX
      if (rank != sizeX) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         rank != sizeX: test failed!" << endl;
        numerr++;
      }
      // MX == M*X
      if (mxcopy != null) {
        tmp = MVT::Clone(*xcopy,sizeX);
        OPT::Apply(*(OM->getOp()),*xcopy,*tmp);
        err = MVDiff(*tmp,*mxcopy);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| MX - M*X ||_F           : " << err << endl;
      }
      // X = xcopy*R
      if (R != null) {
        newR = rcp( new SerialDenseMatrix<int,double>(sizeX,sizeX) );
        OM->innerProd(*xcopy,*X,*newR);
        *newR -= *R;
        err = newR->normFrobenius();
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| newX'*X-R ||_F          : " << err << endl;
      }
      // X^T M X == I
      err = OM->orthonormError(*xcopy);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| X^T M X - I ||_F after  : " << err << endl;
      // X = Q*C + xcopy*R
      if (C != null) {
        tmp = MVT::Clone(*xcopy,sizeX);
        MVT::MvTimesMatAddMv(ONE,*xcopy,*R,ZERO,*tmp);
        MVT::MvTimesMatAddMv(ONE,*Q,*C,ONE,*tmp);
        err = MVDiff(*tmp,*X);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| X-QC-newX ||_F          : " << err << endl;
      }
      // Q^T M X == I
      err = OM->orthogError(*xcopy,*Q);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| Q^T M X ||_F after      : " << err << endl;

      break;
  
    case Undefined:
      sout << "projectAndNormalize() returned Undefined, rank " << rank << endl;

      ind.resize(rank);
      for (int i=0; i<rank; i++) {
        ind[i] = i;
      }
      smlOX = MVT::CloneView(*X,ind);
      smlX = MVT::CloneView(*xcopy,ind); 
      if (mxcopy != null) {
        smlMX = MVT::CloneView(*mxcopy,ind); 
      }
      else {
        smlMX = smlX;
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
        smlR = rcp( new SerialDenseMatrix<int,double>(Teuchos::Copy,*R,rank,rank) );
        newR = rcp( new SerialDenseMatrix<int,double>(rank,rank) );
        OM->innerProd(*smlX,*smlOX,*newR);
        *newR -= *smlR;
        err = newR->normFrobenius();
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| newX'*X-R ||_F          : " << err << endl;
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
        smlR = rcp( new SerialDenseMatrix<int,double>(Teuchos::Copy,*R,rank,rank) );
        smlC = rcp( new SerialDenseMatrix<int,double>(Teuchos::Copy,*C,sizeQ,rank) );
        tmp = MVT::Clone(*smlX,rank);
        MVT::MvTimesMatAddMv(ONE,*smlX,*smlR,ZERO,*tmp);
        MVT::MvTimesMatAddMv(ONE,*Q,*smlC,ONE,*tmp);
        err = MVDiff(*tmp,*smlOX);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| X-QC-newX ||_F          : " << err << endl;
      }
      // Q^T M X == I
      err = OM->orthogError(*smlX,*Q);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| Q^T M X ||_F after      : " << err << endl;
      break;

    case Failed:
      sout << "   -------------------------------------------         projectAndNormalize() Failed" << endl;
      numerr++;
      break;
    default:   
      sout << "   -------------------------------------------         projectAndNormalize() returned invalid value" << endl;
      numerr++;
      break;
    }

  } // test for

  if (verbose || numerr>0) {
    cout << sout.str();
    cout << endl;
  }

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testProject(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                RefCountPtr<const MV> X, RefCountPtr<const MV> Q) {

  const int sizeX = MVT::GetNumberVecs(*X);
  const int sizeQ = MVT::GetNumberVecs(*Q);
  double err;
  ReturnType ret;
  RefCountPtr<MV> tmp;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;
  std::ostringstream sout;
  RefCountPtr<MV> xcopy, mxcopy;
  RefCountPtr<SerialDenseMatrix<int,double> > C;

  int numtests;
  RefCountPtr<MV> MX;
  if (hasM) {
    numtests = 4;
    MX = MVT::Clone(*X,sizeX);
    if ( OPT::Apply(*(OM->getOp()),*X,*MX) != Ok ) return 1;
  }
  else {
    MX = MVT::CloneCopy(*X);
    numtests = 2;
  }

  // test ortho error before orthonormalizing
  err = OM->orthogError(*Q,*X);
  sout << "   || Q^T M X ||_F before     : " << err << endl;

  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
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
      C = rcp( new SerialDenseMatrix<int,double>(sizeQ,sizeX) );
    }
    else {
      C = null;
    }

    ret = OM->project(*xcopy,mxcopy,C,*Q);
    switch (ret) {
    case Ok:
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
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| X-QC-newX ||_F          : " << err << endl;
      }
      // Q^T M X == I
      err = OM->orthogError(*xcopy,*Q);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| Q^T M X ||_F after      : " << err << endl;
      break;
  
    case Failed:
      sout << "   -------------------------------------------         project() Failed" << endl;
      numerr++;
      break;
  
    default:   
      sout << "   -------------------------------------------         project() returned invalid value" << endl;
      numerr++;
      break;
    }

  } // for test

  if (verbose || numerr>0) {
    cout << sout.str();
    cout << endl;
  }

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X) {

  const int sizeX = MVT::GetNumberVecs(*X);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;
  std::ostringstream sout;
  RefCountPtr<MV> xcopy, mxcopy;
  RefCountPtr<SerialDenseMatrix<int,double> > R;

  int numtests;
  RefCountPtr<MV> MX;
  if (hasM) {
    numtests = 4;
    MX = MVT::Clone(*X,sizeX);
    if ( OPT::Apply(*(OM->getOp()),*X,*MX) != Ok ) return 1;
  }
  else {
    MX = MVT::CloneCopy(*X);
    numtests = 2;
  }


  // test ortho error before orthonormalizing
  err = OM->orthonormError(*X);
  sout << "   || X^T M X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
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
      R = rcp( new SerialDenseMatrix<int,double>(sizeX,sizeX) );
      R->putScalar(ZERO);
    }
    else {
      R = null;
    }

    ret = OM->normalize(*xcopy,mxcopy,R,rank);
    switch (ret) {
    case Ok:
      // check the following:
      // rank == orig
      if (rank != sizeX) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         rank != orig: test failed!" << endl;
        numerr++;
      }
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
      // X = xcopy*R
      if (R != null) {
        tmp = MVT::Clone(*xcopy,sizeX);
        MVT::MvTimesMatAddMv(ONE,*xcopy,*R,ZERO,*tmp);
        err = MVDiff(*tmp,*X);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| newX'*X-R ||_F          : " << err << endl;
      }
      // X^T M X == I
      err = OM->orthonormError(*xcopy);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "  " << t << "|| X^T M X - I ||_F after  : " << err << endl;
      break;
  
    case Undefined:
      sout << "normalize() returned Undefined, rank " << rank << endl;
  
      ind.resize(rank);
      for (int i=0; i<rank; i++) {
        ind[i] = i;
      }
      smlX = MVT::CloneView(*xcopy,ind); 
      if (mxcopy != null) {
        smlMX = MVT::CloneView(*mxcopy,ind); 
      }
      else {
        smlMX = smlX;
      }

      // MX == M*X
      if (smlMX != null) {
        tmp = MVT::CloneCopy(*smlX);
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
        SerialDenseMatrix<int,double> smlR(Teuchos::Copy,*R,rank,rank);
        tmp = MVT::Clone(*smlX,rank);
        MVT::MvTimesMatAddMv(ONE,*smlX,smlR,ZERO,*tmp);
        err = MVDiff(*tmp,*X);
        if (err > TOL) {
          sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
          numerr++;
        }
        sout << "  " << t << "|| newX'*X-R ||_F          : " << err << endl;
      }
      // X^T M X == I
      err = OM->orthonormError(*smlX);
      if (err > TOL) {
        sout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      sout << "   || X^T M X - I ||_F after  : " << err << endl;
      break;
    case Failed:
      sout << "   -------------------------------------------         normalize() Failed" << endl;
      numerr++;
      break;
    default:   
      sout << "   -------------------------------------------         normalize() returned invalid value" << endl;
      numerr++;
      break;
    }

  } // for test

  if (verbose || numerr>0) {
    cout << sout.str();
    cout << endl;
  }

  return numerr;
}


double MVDiff(const MV &X, const MV &Y) {
  const int sizeX = MVT::GetNumberVecs(X);
  SerialDenseMatrix<int,double> xTmx(sizeX,sizeX);

  // tmp <- X
  RefCountPtr<MV> tmp = MVT::CloneCopy(X);

  // tmp <- tmp - Y
  MVT::MvAddMv(-ONE,Y,ONE,*tmp,*tmp);

  MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
  double err = 0;
  for (int i=0; i<sizeX; i++) {
    err += xTmx(i,i);
  }
  return SCT::squareroot(err);
}
