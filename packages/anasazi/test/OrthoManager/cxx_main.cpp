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


double orthError(const MV &X, const MV &MX) {
  int rank = MVT::GetNumberVecs(X);
  SerialDenseMatrix<int,double> xTmx(rank,rank);
  MVT::MvTransMv(ONE,X,MX,xTmx);
  for (int i=0; i<rank; i++) {
    xTmx(i,i) -= ONE;
  }
  return xTmx.normFrobenius();
}

int testProject(RefCountPtr<OrthoManager<double,MV,OP> > OM, 
                  MV &X, MV &MX, const OP *M, int orig) {

  SerialDenseMatrix<int,double> xTmx(orig,orig);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;


  // test ortho error before orthonormalizing
  err = orthError(X,MX);
  cout << "   || X^T M X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy = MVT::CloneCopy(X),
                  mxcopy = MVT::CloneCopy(MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,M,*xcopy,orig,2);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  err = orthError(*xcopy,*mxcopy);
  cout << "       MSU.massOrthonormalize error: " << err << endl;


  ret = OM->normalize(X,MX,M,rank);
  switch (ret) {
  case Ok:
    cout << "   normalize() returned Ok" << endl;
    // check the following:
    // rank == orig
    if (rank != orig) {
      cout << "   rank != orig: test failed!" << endl;
      return 1;
    }
    // MX = M*X
    tmp = MVT::CloneCopy(X);
    OPT::Apply(*M,X,*tmp);
    MVT::MvAddMv(-ONE,MX,ONE,*tmp,*tmp);
    MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
    err = 0;
    for (int i=0; i<orig; i++) {
      err += xTmx(i,i);
    }
    err = SCT::squareroot(err);
    cout << "   || MX - M*X ||_F           : " << err << endl;
    // X^T M X == I
    err = orthError(X,MX);
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
    smlX = MVT::CloneView(X,ind); 
    smlMX = MVT::CloneView(MX,ind); 

    xTmx.shape(rank,rank);

    // check the following:
    // MX = M*X
    tmp = MVT::CloneCopy(X);
    OPT::Apply(*M,X,*tmp);
    MVT::MvAddMv(-ONE,MX,ONE,*tmp,*tmp);
    MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
    err = 0;
    for (int i=0; i<orig; i++) {
      err += xTmx(i,i);
    }
    err = SCT::squareroot(err);
    cout << "   || MX - M*X ||_F           : " << err << endl;
    // X^T M X == I
    err = orthError(*smlX,*smlMX);
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Failed:
    cout << "normalize() Failed" << endl;
    return 1;

  default:   
    cout << "normalize() returned invalid value" << endl;
    return 1;
  }

  return 0;
}




int testNormalize(RefCountPtr<OrthoManager<double,MV,OP> > OM, 
                  MV &X, MV &MX, const OP *M, int orig) {

  SerialDenseMatrix<int,double> xTmx(orig,orig);
  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> tmp, smlX, smlMX;
  std::vector<int> ind;


  // test ortho error before orthonormalizing
  err = orthError(X,MX);
  cout << "   || X^T M X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy = MVT::CloneCopy(X),
                  mxcopy = MVT::CloneCopy(MX);
  int iret = MSU.massOrthonormalize(*xcopy,*mxcopy,M,*xcopy,orig,2);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  err = orthError(*xcopy,*mxcopy);
  cout << "       MSU.massOrthonormalize error: " << err << endl;


  ret = OM->normalize(X,MX,M,rank);
  switch (ret) {
  case Ok:
    cout << "   normalize() returned Ok" << endl;
    // check the following:
    // rank == orig
    if (rank != orig) {
      cout << "   rank != orig: test failed!" << endl;
      return 1;
    }
    // MX = M*X
    tmp = MVT::CloneCopy(X);
    OPT::Apply(*M,X,*tmp);
    MVT::MvAddMv(-ONE,MX,ONE,*tmp,*tmp);
    MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
    err = 0;
    for (int i=0; i<orig; i++) {
      err += xTmx(i,i);
    }
    err = SCT::squareroot(err);
    cout << "   || MX - M*X ||_F           : " << err << endl;
    // X^T M X == I
    err = orthError(X,MX);
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
    smlX = MVT::CloneView(X,ind); 
    smlMX = MVT::CloneView(MX,ind); 

    xTmx.shape(rank,rank);

    // check the following:
    // MX = M*X
    tmp = MVT::CloneCopy(X);
    OPT::Apply(*M,X,*tmp);
    MVT::MvAddMv(-ONE,MX,ONE,*tmp,*tmp);
    MVT::MvTransMv(ONE,*tmp,*tmp,xTmx);
    err = 0;
    for (int i=0; i<orig; i++) {
      err += xTmx(i,i);
    }
    err = SCT::squareroot(err);
    cout << "   || MX - M*X ||_F           : " << err << endl;
    // X^T M X == I
    err = orthError(*smlX,*smlMX);
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Failed:
    cout << "normalize() Failed" << endl;
    return 1;

  default:   
    cout << "normalize() returned invalid value" << endl;
    return 1;
  }

  return 0;
}





int testNormalize(RefCountPtr<OrthoManager<double,MV,OP> > OM, 
                  MV &X, int orig) {

  double err;
  int rank;
  ReturnType ret;
  RefCountPtr<MV> smlX;
  std::vector<int> ind;


  // test ortho error before orthonormalizing
  err = orthError(X,X);
  cout << "   || X^T X - I ||_F before : " << err << endl;


  // first, run with MSUtils
  ModalSolverUtils<double,MV,OP> MSU(rcp(new OutputManager<double>()));
  RefCountPtr<MV> xcopy = MVT::CloneCopy(X);
  int iret = MSU.massOrthonormalize(*xcopy,*xcopy,NULL,*xcopy,orig,2);
  cout << "       MSU.massOrthonormalize returned " << iret << endl;
  err = orthError(*xcopy,*xcopy);
  cout << "       MSU.massOrthonormalize error: " << err << endl;




  ret = OM->normalize(X,X,0,rank);
  switch (ret) {
  case Ok:
    cout << "   normalize() returned Ok" << endl;
    // check the following:
    // rank == orig
    if (rank != orig) {
      cout << "   rank != orig: test failed!" << endl;
      return 1;
    }
    // X^T X == I
    err = orthError(X,X);
    cout << "   || X^T X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Undefined:
    cout << "normalize() returned Undefined, rank " << rank << endl;

    ind.resize(rank);
    for (int i=0; i<rank; i++) {
      ind[i] = i;
    }
    smlX = MVT::CloneView(X,ind); 

    // check the following:
    // X^T X == I
    err = orthError(*smlX,*smlX);
    cout << "   || X^T M X - I ||_F after  : " << err << endl;
    // span(X before) \subspace span(X after)
    // FINISH: test canonical angles between oldX and newX ?

    break;

  case Failed:
    cout << "normalize() Failed" << endl;
    return 1;

  default:   
    cout << "normalize() returned invalid value" << endl;
    return 1;
  }

  return 0;
}


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
    cout << Anasazi_Version() << endl << endl;
  }
  
  ReturnType ret;
  
  // Problem information
  const int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;
  const int orig = 10;

  // Create problem
  RefCountPtr<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  // Get the mass matrix
  RefCountPtr<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );

  // Create a basic ortho manager 
  RefCountPtr<OrthoManager<double,MV,OP> > BOM = rcp( new BasicOrthoManager<double,MV,OP>() );

  // multivector to spawn off of
  RefCountPtr<MV> X = rcp( new Epetra_MultiVector(M->OperatorDomainMap(), orig) );
  RefCountPtr<MV> MX = MVT::CloneCopy(*X);

  {
    // check normalize() on a random multivector with M inner product
    cout << endl;
    cout << " normalize() with M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    OPT::Apply(*M,*X,*MX);
    numFailed += testNormalize(BOM,*X,*MX,M.get(),orig);
  }
  {
    // check normalize() on a random multivector with Euclidean inner product
    cout << endl;
    cout << " normalize() without M testing on random multivector " << endl;
    MVT::MvRandom(*X);
    numFailed += testNormalize(BOM,*X,orig);
  }


  {
    // check normalize() on a rank-deficient multivector with M inner product
    cout << endl;
    cout << " normalize() with M testing on rank-deficient multivector " << endl;
    std::vector<int> ind(1); 
    // Assuming here that orig == 10
    MVT::MvRandom(*X);
    ind[0] = 2;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    ind[0] = 7;
    MVT::SetBlock(*mid,ind,*X);
    
    OPT::Apply(*M,*X,*MX);
    
    numFailed += testNormalize(BOM,*X,*MX,M.get(),orig);
  }
  {
    // check normalize() on a rank-deficient multivector with Euclidean inner product
    cout << endl;
    cout << " normalize() without M testing on rank-deficient multivector " << endl;
    std::vector<int> ind(1); 
    // Assuming here that orig == 10
    MVT::MvRandom(*X);
    ind[0] = 2;
    RefCountPtr<MV> mid = MVT::CloneCopy(*X,ind);
    ind[0] = 7;
    MVT::SetBlock(*mid,ind,*X);
    
    numFailed += testNormalize(BOM,*X,orig);
  }


  cout << endl;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (numFailed) {
    if (verbose && MyPID==0)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;
}	
