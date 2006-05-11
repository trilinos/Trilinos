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
int testProjectAndNormalize(RefCountPtr<OrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X, 
                            RefCountPtr<MV> MX, RefCountPtr<const OP> M, RefCountPtr<const MV> Q);
int testProject(RefCountPtr<OrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X, 
                RefCountPtr<MV> MX, RefCountPtr<const OP> M, RefCountPtr<const MV> Q);
int testNormalize(RefCountPtr<OrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X, 
                  RefCountPtr<MV> MX, RefCountPtr<const OP> M);
double orthonormError(const MV &X1, const MV &X2);
double orthogError(const MV &X1, const MV &X2);


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
  RefCountPtr<OrthoManager<double,MV,OP> > BOM = rcp( new BasicOrthoManager<double,MV,OP>() );

  // multivector to spawn off of
  RefCountPtr<MV> X = rcp( new Epetra_MultiVector(M->OperatorDomainMap(), sizeX) );
  RefCountPtr<MV> MX = MVT::CloneCopy(*X);

  cout << endl;
  cout << " Generating Q1,MQ1 for project() testing... " << endl;
  RefCountPtr<MV> Q1 = MVT::Clone(*X,sizeQ),
                 MQ1 = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q1);
  OPT::Apply(*M,*Q1,*MQ1);
  int dummy = 0;
  if ( BOM->normalize(*Q1,*MQ1,M.get(),dummy) != Ok ) {
    cout << " FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
  }
  err = orthonormError(*Q1,*MQ1);
  if (err > TOL) {
    cout << " FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
  }
  cout << "   || Q1^T M Q1 - I ||_F        : " << err << endl;

  cout << " Generating Q2 for project() testing... " << endl;
  RefCountPtr<MV> Q2 = MVT::Clone(*X,sizeQ);
  MVT::MvRandom(*Q2);
  if ( BOM->normalize(*Q2,*Q2,NULL,dummy) != Ok ) {
    cout << " FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
  }
  err = orthonormError(*Q2,*Q2);
  if (err > TOL) {
    cout << " FAILED!!! Further tests will not be valid." << endl;
    numFailed += 1;
  }
  cout << "   || Q2^T Q2 - I ||_F          : " << err << endl;

  {
    cout << endl;
    cout << " project() with M testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    OPT::Apply(*M,*X,*MX);
    numFailed += testProject(BOM,X,MX,M,Q1);
  }
  {
    cout << endl;
    cout << " project() without M testing on random multivector " << endl;
    
    MVT::MvRandom(*X);
    numFailed += testProject(BOM,X,null,null,Q2);
  }

  {
    cout << endl;
    cout << " normalize() with M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    OPT::Apply(*M,*X,*MX);
    numFailed += testNormalize(BOM,X,MX,M);
  }
  {
    cout << endl;
    cout << " normalize() without M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(BOM,X,null,null);
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
    
    OPT::Apply(*M,*X,*MX);
    
    numFailed += testNormalize(BOM,X,MX,M);
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
    
    numFailed += testNormalize(BOM,X,null,null);
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
    
    OPT::Apply(*M,*X,*MX);
    
    numFailed += testNormalize(BOM,X,MX,M);
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
    
    numFailed += testNormalize(BOM,X,null,null);
  }


  {
    cout << endl;
    cout << " projectAndNormalize() with M testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    OPT::Apply(*M,*X,*MX);
    
    numFailed += testProjectAndNormalize(BOM,X,MX,M,Q1);
  }
  {
    cout << endl;
    cout << " projectAndNormalize() without M testing on random multivector " << endl;
    std::vector<int> ind(1); 
    MVT::MvRandom(*X);
    
    numFailed += testProjectAndNormalize(BOM,X,null,null,Q2);
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
    
    OPT::Apply(*M,*X,*MX);
    
    numFailed += testProjectAndNormalize(BOM,X,MX,M,Q1);
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
    
    numFailed += testProjectAndNormalize(BOM,X,null,null,Q2);
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
    
    OPT::Apply(*M,*X,*MX);
    
    numFailed += testProjectAndNormalize(BOM,X,MX,M,Q1);
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
    
    numFailed += testProjectAndNormalize(BOM,X,null,null,Q1);
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
int testProjectAndNormalize(RefCountPtr<OrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X, 
                            RefCountPtr<MV> MX, RefCountPtr<const OP> M, RefCountPtr<const MV> Q) {

  const int orig = MVT::GetNumberVecs(*X);
  SerialDenseMatrix<int,double> xTmx(orig,orig);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (M != null);
  int numerr = 0;

  if (!hasM) {
    MX = X;
  }

  // test ortho error before orthonormalizing
  err = orthonormError(*X,*MX);
  cout << "   || X^T M X - I ||_F before : " << err << endl;
  err = orthogError(*Q,*MX);
  cout << "   || Q^T M X ||_F before     : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy, mxcopy;
  xcopy = MVT::CloneCopy(*X);
  if (hasM) {
    mxcopy = MVT::CloneCopy(*MX);
  }
  else {
    mxcopy = xcopy;
  }
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,M.get(),*Q,orig,0);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  cout << "       MSU.massOrthonormalize || X^T M X - I ||_F : " << orthonormError(*xcopy,*mxcopy) << endl;
  cout << "       MSU.massOrthonormalize || Q^T M X ||_F     : " << orthogError(*Q,*mxcopy) << endl;

  ret = OM->projectAndNormalize(*X,*MX,M.get(),*Q,rank);
  switch (ret) {
  case Ok:
    cout << "   projectAndNormalize() returned Ok" << endl;
    // check the following:
    // rank == orig
    if (rank != orig) {
      cout << " ----------------------------------------------------- rank != orig: test failed!" << endl;
      numerr++;
    }
    if (hasM) {
    // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*M,*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<orig; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = orthonormError(*X,*MX);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // Q^T M X == I
    err = orthogError(*MX,*Q);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
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
      OPT::Apply(*M,*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<orig; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = orthonormError(*smlX,*smlMX);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?
    // Q^T M X == I
    err = orthogError(*smlMX,*Q);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || Q^T M X ||_F after      : " << err << endl;

    break;

  case Failed:
    cout << " ----------------------------------------------------- projectAndNormalize() Failed" << endl;
    numerr++;

  default:   
    cout << " ----------------------------------------------------- projectAndNormalize() returned invalid value" << endl;
    numerr++;
  }

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testProject(RefCountPtr<OrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X, 
                RefCountPtr<MV> MX, RefCountPtr<const OP> M, RefCountPtr<const MV> Q) {

  int sizeX = MVT::GetNumberVecs(*X);
  SerialDenseMatrix<int,double> xTmx(sizeX,sizeX);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (M != null);
  int numerr = 0;
  
  if (!hasM) {
    MX = X;
  }

  // test ortho error before orthonormalizing
  err = orthogError(*Q,*MX);
  cout << "   || Q^T M X ||_F before     : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy = MVT::CloneCopy(*X),
                  mxcopy = MVT::CloneCopy(*MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,M.get(),*Q,sizeX,1);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  err = orthogError(*Q,*mxcopy);
  cout << "       MSU.massOrthonormalize error: " << err << endl;


  ret = OM->project(*X,*MX,M.get(),*Q);
  switch (ret) {
  case Ok:
    cout << "   project() returned Ok" << endl;
    // check the following:
    if (hasM) {
      // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*M,*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<sizeX; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // Q^T M X == I
    err = orthogError(*MX,*Q);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || Q^T M X ||_F after      : " << err << endl;

    break;

  case Failed:
    cout << " ----------------------------------------------------- project() Failed" << endl;
    numerr++;

  default:   
    cout << " ----------------------------------------------------- project() returned invalid value" << endl;
    numerr++;
  }

  return numerr;
}



//////////////////////////////////////////////////////////////////////
int testNormalize(RefCountPtr<OrthoManager<double,MV,OP> > OM, RefCountPtr<MV> X, 
                  RefCountPtr<MV> MX, RefCountPtr<const OP> M) {

  const int orig = MVT::GetNumberVecs(*X);
  SerialDenseMatrix<int,double> xTmx(orig,orig);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;
  bool hasM = (M != null);
  int numerr = 0;

  if (!hasM) {
    MX = X;
  }

  // test ortho error before orthonormalizing
  err = orthonormError(*X,*MX);
  cout << "   || X^T M X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy, mxcopy;
  xcopy = MVT::CloneCopy(*X);
  if (hasM) {
    mxcopy = MVT::CloneCopy(*MX);
  }
  else {
    mxcopy = xcopy;
  }
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,M.get(),*xcopy,orig,2);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  cout << "       MSU.massOrthonormalize error: " << orthonormError(*xcopy,*mxcopy) << endl;

  ret = OM->normalize(*X,*MX,M.get(),rank);
  switch (ret) {
  case Ok:
    cout << "   normalize() returned Ok" << endl;
    // check the following:
    // rank == orig
    if (rank != orig) {
      cout << " ----------------------------------------------------- rank != orig: test failed!" << endl;
      numerr++;
    }
    if (hasM) {
    // MX = M*X
      tmp = MVT::CloneCopy(*X);
      OPT::Apply(*M,*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<orig; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = orthonormError(*X,*MX);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
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
      OPT::Apply(*M,*X,*tmp);
      MVT::MvAddMv(-ONE,*MX,ONE,*tmp,*tmp);
      MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
      err = 0;
      for (int i=0; i<orig; i++) {
        err += xTmx(i,i);
      }
      err = SCT::squareroot(err);
      if (err > TOL) {
        cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
        numerr++;
      }
      cout << "   || MX - M*X ||_F           : " << err << endl;
    }
    // X^T M X == I
    err = orthonormError(*smlX,*smlMX);
    if (err > TOL) {
      cout << " ----------------------------------------------------- tolerance exceeded! test failed!" << endl;
      numerr++;
    }
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Failed:
    cout << " ----------------------------------------------------- normalize() Failed" << endl;
    numerr++;

  default:   
    cout << " ----------------------------------------------------- normalize() returned invalid value" << endl;
    numerr++;
  }

  return numerr;
}


// compute norm of X1^T X2 - I
double orthonormError(const MV &X1, const MV &X2) {
  int rank = MVT::GetNumberVecs(X1);
  SerialDenseMatrix<int,double> xTx(rank,rank);
  MVT::MvTransMv(ONE,X1,X2,xTx);
  for (int i=0; i<rank; i++) {
    xTx(i,i) -= ONE;
  }
  return xTx.normFrobenius();
}


// compute norm of X1^T X2
double orthogError(const MV &X1, const MV &X2) {
  int r1 = MVT::GetNumberVecs(X1);
  int r2  = MVT::GetNumberVecs(X2);
  SerialDenseMatrix<int,double> xTx(r1,r2);
  MVT::MvTransMv(ONE,X1,X2,xTx);
  return xTx.normFrobenius();
}





