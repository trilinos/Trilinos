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


// some forward declarations
int testProjectAndNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                            RefCountPtr<MV> X, RefCountPtr<const MV> Q);
int testProject(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                RefCountPtr<MV> X, RefCountPtr<const MV> Q);
int testNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X);

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
  bool verbose = true;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  if (verbose && MyPID == 0) {
    cout << Anasazi_Version() << endl;
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

  // Create problem
  RefCountPtr<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  // Get the mass matrix
  RefCountPtr<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  // Create a basic ortho manager 
  RefCountPtr<MatOrthoManager<double,MV,OP> > OM   = rcp( new BasicOrthoManager<double,MV,OP>() );
  RefCountPtr<MatOrthoManager<double,MV,OP> > OM_M = rcp( new BasicOrthoManager<double,MV,OP>(M) );

  // multivector to spawn off of
  RefCountPtr<MV> X = rcp( new Epetra_MultiVector(M->OperatorDomainMap(), sizeX) );

  cout << endl;
  cout << " Generating Q1 for project() testing... " << endl;
  RefCountPtr<MV> Q1  = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q1);
  int dummy = 0;
  ret = OM_M->normalize(*Q1,null,dummy);
  if ( ret  != Ok ) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! normalize() return " << ret << "/" << dummy << " . Further tests will not be valid." << endl;
    numFailed += 1;
  }
  err = OM_M->orthonormError(*Q1);
  if (err > TOL) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
  }
  cout << "   || Q1^T M Q1 - I ||_F        : " << err << endl;

  cout << " Generating Q2 for project() testing... " << endl;
  RefCountPtr<MV> Q2 = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q2);
  ret = OM->normalize(*Q2,null,dummy);
  if ( ret != Ok ) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! normalize() return " << ret << "/" << dummy << " . Further tests will not be valid." << endl;
    numFailed += 1;
  }
  err = OM->orthonormError(*Q2);
  if (err > TOL) {
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
  }
  cout << "   || Q2^T Q2 - I ||_F          : " << err << endl;

  {
    cout << endl;
    cout << " project() with M testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(OM_M,X,Q1);
  }
  {
    cout << endl;
    cout << " project() without M testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(OM,X,Q2);
  }

  {
    cout << endl;
    cout << " normalize() with M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(OM_M,X);
  }
  {
    cout << endl;
    cout << " normalize() without M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(OM,X);
  }


  {
    cout << endl;
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
    cout << endl;
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
    cout << endl;
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
    cout << endl;
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
    cout << endl;
    cout << " projectAndNormalize() with M testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(OM_M,X,Q1);
  }
  {
    cout << endl;
    cout << " projectAndNormalize() without M testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(OM,X,Q2);
  }

  {
    cout << endl;
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
    cout << endl;
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
    cout << endl;
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
    cout << endl;
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

  cout << endl;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (numFailed) {
    if (verbose && MyPID==0) {
      cout << numFailed << " errors." << endl;
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;
}	




//////////////////////////////////////////////////////////////////////
int testProjectAndNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                            RefCountPtr<MV> X, RefCountPtr<const MV> Q) {

  const int sizeX = MVT::GetNumberVecs(*X);
  SerialDenseMatrix<int,double> xTmx(sizeX,sizeX);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX;
  std::vector<int> ind;
  int numerr = 0;
  bool hasM = (OM->getOp() != null);


  // test ortho error before orthonormalizing
  err = OM->orthonormError(*X);
  cout << "   || X^T M X - I ||_F before : " << err << endl;
  err = OM->orthogError(*X,*Q);
  cout << "   || Q^T M X ||_F before     : " << err << endl;

  RefCountPtr<MV> MX = MVT::Clone(*X,sizeX);
  if (hasM) {
    if ( OPT::Apply(*(OM->getOp()),*X,*MX) != Ok ) return 1;
  }
  // else, it won't be referenced

  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy = MVT::CloneCopy(*X),
                 mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,OM->getOp().get(),*Q,sizeX,0);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  cout << "       MSU.massOrthonormalize || X^T M X - I ||_F : " << OM->orthonormError(*xcopy) << endl;
  cout << "       MSU.massOrthonormalize || Q^T M X ||_F     : " << OM->orthogError(*Q,*xcopy) << endl;

  ret = OM->projectAndNormalize(*X,MX,null,null,*Q,rank);
  switch (ret) {
  case Ok:
    cout << "   projectAndNormalize() returned Ok" << endl;
    // check the following:
    // rank == sizeX
    if (rank != sizeX) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         rank != sizeX: test failed!" << endl;
      numerr++;
    }
    if (hasM) {
    // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*(OM->getOp()),*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<sizeX; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    // FINISH: test both orthonormError methods
    err = OM->orthonormError(*X);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // Q^T M X == I
    // FINISH: test both orthogError methods
    err = OM->orthogError(*X,*Q);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || Q^T M X ||_F after      : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Undefined:
    cout << "projectAndNormalize() returned Undefined, rank " << rank << endl;

    ind.resize(rank);
    for (int i=0; i<rank; i++) {
      ind[i] = i;
    }
    smlX = MVT::CloneView(*X,ind); 
    xTmx.shape(rank,rank);

    // check the following:
    if (OM->getOp() != Teuchos::null) {
      // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*(OM->getOp()),*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<sizeX; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = OM->orthonormError(*smlX);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?
    // Q^T M X == I
    err = OM->orthogError(*smlX,*Q);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || Q^T M X ||_F after      : " << err << endl;

    break;

  case Failed:
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         projectAndNormalize() Failed" << endl;
    numerr++;

  default:   
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         projectAndNormalize() returned invalid value" << endl;
    numerr++;
  }

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testProject(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, 
                RefCountPtr<MV> X, RefCountPtr<const MV> Q) {

  const int sizeX = MVT::GetNumberVecs(*X);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;
  
  // test ortho error before orthonormalizing
  err = OM->orthogError(*Q,*X);
  cout << "   || Q^T M X ||_F before     : " << err << endl;

  RefCountPtr<MV> MX = MVT::Clone(*X,sizeX);
  if (hasM) {
    if ( OPT::Apply(*(OM->getOp()),*X,*MX) != Ok ) return 1;
  }
  // else, it won't be referenced

  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy = MVT::CloneCopy(*X),
                  mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,OM->getOp().get(),*Q,sizeX,1);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  err = OM->orthogError(*Q,*xcopy);
  cout << "       MSU.massOrthonormalize error: " << err << endl;


  ret = OM->project(*X,MX,null,*Q);
  switch (ret) {
  case Ok:
    cout << "   project() returned Ok" << endl;
    // check the following:
    if (hasM) {
      // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*(OM->getOp()),*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      SerialDenseMatrix<int,double> xTmx(sizeX,sizeX);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<sizeX; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // Q^T M X == I
    err = OM->orthogError(*X,*Q);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || Q^T M X ||_F after      : " << err << endl;

    break;

  case Failed:
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         project() Failed" << endl;
    numerr++;

  default:   
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         project() returned invalid value" << endl;
    numerr++;
  }

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testNormalize(RefCountPtr<MatOrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X) {

  const int sizeX = MVT::GetNumberVecs(*X);
  SerialDenseMatrix<int,double> xTmx(sizeX,sizeX);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (OM->getOp() != null);
  int numerr = 0;

  RefCountPtr<MV> MX = MVT::Clone(*X,sizeX);
  if (hasM) {
    if ( OPT::Apply(*(OM->getOp()),*X,*MX) != Ok ) return 1;
  }
  // else, it won't be referenced


  // test ortho error before orthonormalizing
  err = OM->orthonormError(*X);
  cout << "   || X^T M X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy  = MVT::CloneCopy(*X),
                  mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,OM->getOp().get(),*xcopy,sizeX,2);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  cout << "       MSU.massOrthonormalize error: " << OM->orthonormError(*xcopy) << endl;

  ret = OM->normalize(*X,MX,null,rank);
  switch (ret) {
  case Ok:
    cout << "   normalize() returned Ok" << endl;
    // check the following:
    // rank == orig
    if (rank != sizeX) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         rank != orig: test failed!" << endl;
      numerr++;
    }
    if (hasM) {
    // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*(OM->getOp()),*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<sizeX; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = OM->orthonormError(*X);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Undefined:
    cout << "normalize() returned Undefined, rank " << rank << endl;

    ind.resize(rank);
    for (int i=0; i<rank; i++) {
      ind[i] = i;
    }
    smlX = MVT::CloneView(*X,ind); 
    if (hasM) {
      smlMX = MVT::CloneView(*MX,ind); 
    }
    else {
      smlMX = smlX;
    }
    xTmx.shape(rank,rank);

    // check the following:
    if (hasM) {
      // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*(OM->getOp()),*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<sizeX; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = OM->orthonormError(*smlX);
    if (err > TOL) {
      cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Failed:
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         normalize() Failed" << endl;
    numerr++;

  default:   
    cout << "   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv         normalize() returned invalid value" << endl;
    numerr++;
  }

  return numerr;
}

