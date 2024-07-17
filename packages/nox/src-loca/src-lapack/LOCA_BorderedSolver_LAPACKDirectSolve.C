// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_LAPACKDirectSolve.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_LAPACK.H"
#include "Teuchos_LAPACK.hpp"
#include "LOCA_BorderedSolver_JacobianOperator.H"
#include "LOCA_BorderedSolver_ComplexOperator.H"
#include "LOCA_Hopf_ComplexMultiVector.H"

LOCA::BorderedSolver::LAPACKDirectSolve::LAPACKDirectSolve(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& slvrParams):
  globalData(global_data),
  solverParams(slvrParams),
  grp(),
  op(),
  A(),
  B(),
  C(),
  augJacSolver(),
#ifdef HAVE_TEUCHOS_COMPLEX
  augComplexSolver(),
#endif
  n(0),
  m(0),
  N(0),
  isZeroA(true),
  isZeroB(true),
  isZeroC(true),
  isZeroF(true),
  isZeroG(true),
  isComplex(false)
{
}

LOCA::BorderedSolver::LAPACKDirectSolve::~LAPACKDirectSolve()
{
}

  void
LOCA::BorderedSolver::LAPACKDirectSolve::setMatrixBlocks(
    const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& oper,
    const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
    const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
    const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC)
{
  std::string callingFunction =
    "LOCA::BorderedSolver::LAPACKDirectSolve::setMatrixBlocks()";

  // Set block pointers
  op = oper;
  A = blockA;

  B = Teuchos::rcp_dynamic_cast<const LOCA::MultiContinuation::ConstraintInterfaceMVDX>(blockB);
  if (B.get() == NULL)
    globalData->locaErrorCheck->throwError(
        callingFunction,
        std::string("Constraint argument is not of type\n") +
        std::string("LOCA::MultiContinuation::ConstraintInterfaceMVDX!\n") +
        std::string("The LAPACK Direct Solve bordered solver method can\n") +
        std::string("only be used with constraints that support obtaining\n") +
        std::string("the constraint derivative as a multivector."));

  C = blockC;

  isZeroA = (A.get() == NULL);
  isZeroB = B->isDXZero();
  isZeroC = (C.get() == NULL);

  // ensure blocks B and C are not both zero
  if (isZeroB && isZeroC)
    globalData->locaErrorCheck->throwError(
        callingFunction,
        "Blocks B and C cannot both be zero");

  // ensure blocks A and C are not both zero
  if (isZeroA && isZeroC)
    globalData->locaErrorCheck->throwError(
        callingFunction,
        "Blocks A and C cannot both be zero");

  // Fill augmented matrix
  Teuchos::RCP<const LOCA::BorderedSolver::JacobianOperator> jacOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::JacobianOperator>(op);
  Teuchos::RCP<const LOCA::BorderedSolver::ComplexOperator> complexOp =
    Teuchos::rcp_dynamic_cast<const LOCA::BorderedSolver::ComplexOperator>(op);

  if (jacOp != Teuchos::null) {
    const NOX::LAPACK::Vector *v;
    const NOX::Abstract::MultiVector *BV;

    isComplex = false;

    grp = Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(jacOp->getGroup());
    if (grp.get() == NULL)
      globalData->locaErrorCheck->throwError(
          callingFunction,
          std::string("Group argument is not of type LOCA::LAPACK::Group!\n") +
          std::string("The LAPACK Direct Solve bordered solver method can\n") +
          std::string("only be used with LAPACK groups."));

    // Get the Jacobian matrix and size
    const NOX::LAPACK::Matrix<double>& J = grp->getJacobianMatrix();
    n = J.numRows();

    // Get the number of additional rows/columns
    if (!isZeroA)
      m = A->numVectors();
    else
      m = C->numCols();

    // Form a new (n+m) x (n+m) matrix if this is a new size
    if (n+m != N) {
      N = n+m;
      augJacSolver = Teuchos::rcp(new NOX::LAPACK::LinearSolver<double>(N));
    }
    else {
      augJacSolver->reset();
    }
    NOX::LAPACK::Matrix<double>& augmentedJ = augJacSolver->getMatrix();

    // Copy Jacobian
    for (int j=0; j<n; j++)
      for (int i=0; i<n; i++)
        augmentedJ(i,j) = J(i,j);

    // Copy A
    if (isZeroA) {
      for (int j=0; j<m; j++)
        for (int i=0; i<n; i++)
          augmentedJ(i,j+n) = 0.0;
    }
    else {
      for (int j=0; j<m; j++) {
        v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*A)[j]);
        TEUCHOS_ASSERT(v != NULL);
        for (int i=0; i<n; i++)
          augmentedJ(i,j+n) = (*v)(i);
      }
    }

    // Copy B
    if (isZeroB) {
      for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
          augmentedJ(i+n,j) = 0.0;
    }
    else {
      BV = B->getDX();
      for (int i=0; i<m; i++) {
        v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*BV)[i]);
        TEUCHOS_ASSERT(v != NULL);
        for (int j=0; j<n; j++)
          augmentedJ(i+n,j) = (*v)(j);
      }
    }

    // Copy C
    if (isZeroC) {
      for (int j=0; j<m; j++)
        for (int i=0; i<m; i++)
          augmentedJ(i+n,j+n) = 0.0;
    }
    else {
      for (int j=0; j<m; j++)
        for (int i=0; i<m; i++)
          augmentedJ(i+n,j+n) = (*C)(i,j);
    }

  }
  else if (complexOp != Teuchos::null) {
#ifdef HAVE_TEUCHOS_COMPLEX
    Teuchos::RCP<const LOCA::Hopf::ComplexMultiVector> cA;
    Teuchos::RCP<const LOCA::Hopf::ComplexMultiVector> cB;
    Teuchos::RCP<const NOX::Abstract::MultiVector> A_real;
    Teuchos::RCP<const NOX::Abstract::MultiVector> A_imag;
    Teuchos::RCP<const NOX::Abstract::MultiVector> B_real;
    Teuchos::RCP<const NOX::Abstract::MultiVector> B_imag;
    const NOX::LAPACK::Vector *v1, *v2;
    Teuchos::RCP<const NOX::Abstract::MultiVector> BV;

    isComplex = true;

    if (!isZeroA) {
      cA =
        Teuchos::rcp_dynamic_cast<const LOCA::Hopf::ComplexMultiVector>(A,
            true);
      A_real = cA->getRealMultiVec();
      A_imag = cA->getImagMultiVec();
    }
    if (!isZeroB) {
      BV = Teuchos::rcp(B->getDX(),false);
      cB =
        Teuchos::rcp_dynamic_cast<const LOCA::Hopf::ComplexMultiVector>(BV,
            true);
      B_real = cB->getRealMultiVec();
      B_imag = cB->getImagMultiVec();
    }

    grp = Teuchos::rcp_dynamic_cast<const LOCA::LAPACK::Group>(complexOp->getGroup());
    if (grp.get() == NULL)
      globalData->locaErrorCheck->throwError(
          callingFunction,
          std::string("Group argument is not of type LOCA::LAPACK::Group!\n") +
          std::string("The LAPACK Direct Solve bordered solver method can\n") +
          std::string("only be used with LAPACK groups."));

    // Get the number of additional rows/columns
    if (!isZeroA)
      m = A->numVectors()/2;  // Two columns for each complex vector
    else
      m = C->numCols()/2;

    // Get the complex matrix and size
    const NOX::LAPACK::Matrix< std::complex<double> >& mat =
      grp->getComplexMatrix();
    n = mat.numRows();

    // Form a new (n+m) x (n+m) matrix if this is a new size
    if (n+m != N) {
      N = n+m;
      augComplexSolver =
        Teuchos::rcp(new NOX::LAPACK::LinearSolver< std::complex<double> >(N));
    }
    else {
      augComplexSolver->reset();
    }
    NOX::LAPACK::Matrix< std::complex<double> >& augmentedMat =
      augComplexSolver->getMatrix();

    // Copy matrix
    for (int j=0; j<n; j++)
      for (int i=0; i<n; i++)
        augmentedMat(i,j) = mat(i,j);

    // Copy A
    if (isZeroA) {
      for (int j=0; j<m; j++)
        for (int i=0; i<n; i++)
          augmentedMat(i,j+n) = 0.0;
    }
    else {
      for (int j=0; j<m; j++) {
        v1 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*A_real)[2*j]);
        TEUCHOS_ASSERT(v1 != NULL);
        v2 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*A_imag)[2*j]);
        TEUCHOS_ASSERT(v2 != NULL);
        for (int i=0; i<n; i++)
          augmentedMat(i,j+n) = std::complex<double>((*v1)(i), (*v2)(i));
      }
    }

    // Copy B
    if (isZeroB) {
      for (int i=0; i<m; i++)
        for (int j=0; j<n; j++)
          augmentedMat(i+n,j) = 0.0;
    }
    else {
      for (int i=0; i<m; i++) {
        v1 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*B_real)[2*i]);
        TEUCHOS_ASSERT(v1 != NULL);
        v2 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*B_imag)[2*i]);
        TEUCHOS_ASSERT(v2 != NULL);
        for (int j=0; j<n; j++)
          augmentedMat(i+n,j) = std::complex<double>((*v1)(j), -(*v2)(j));
      }
    }

    // Copy C
    if (isZeroC) {
      for (int j=0; j<m; j++)
        for (int i=0; i<m; i++)
          augmentedMat(i+n,j+n) = 0.0;
    }
    else {
      for (int j=0; j<m; j++)
        for (int i=0; i<m; i++)
          augmentedMat(i+n,j+n) = std::complex<double>((*C)(i,2*j),
              (*C)(i+m,2*j));
    }
#else
    globalData->locaErrorCheck->throwError(
        callingFunction,
        "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
#endif
  }

  else {
    globalData->locaErrorCheck->throwError(
        callingFunction,
        std::string("Op argument must be of type !\n") +
        std::string("LOCA::BorderedSolver::JacobianOperator or \n") +
        std::string("LOCA::BorderedSolver::ComplexOperator."));
  }
}

  NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::initForSolve()
{
  return NOX::Abstract::Group::Ok;
}

  NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::initForTransposeSolve()
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::apply(
    const NOX::Abstract::MultiVector& X,
    const NOX::Abstract::MultiVector::DenseMatrix& Y,
    NOX::Abstract::MultiVector& U,
    NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  //   int numCols = X.numVectors();

  //   if (!isComplex) {
  //     const NOX::LAPACK::Vector *v;
  //     NOX::LAPACK::Vector *w;

  //     // Concatenate X & Y into a single matrix
  //     NOX::LAPACK::Matrix<double> RHS(N,numCols);
  //     NOX::LAPACK::Matrix<double> LHS(N,numCols);
  //     for (int j=0; j<numCols; j++) {
  //       v = dynamic_cast<const NOX::LAPACK::Vector*>(&X[j]);
  //       for (int i=0; i<n; i++)
  //     RHS(i,j) = (*v)(i);
  //       for (int i=0; i<m; i++)
  //     RHS(i+n,j) = Y(i,j);
  //     }

  //     // Solve for LHS
  //     augJacSolver->apply(false, numCols, &RHS(0,0), &LHS(0,0));

  //     // Copy result into U and V
  //     for (int j=0; j<numCols; j++) {
  //       w = dynamic_cast<NOX::LAPACK::Vector*>(&U[j]);
  //       for (int i=0; i<n; i++)
  //     (*w)(i) = LHS(i,j);
  //       for (int i=0; i<m; i++)
  //     V(i,j) = LHS(n+i,j);
  //     }
  //   }

  //   else {
  //     const LOCA::Hopf::ComplexMultiVector* cX;
  //     LOCA::Hopf::ComplexMultiVector* cU;
  //     Teuchos::RCP<const NOX::Abstract::MultiVector> X_real;
  //     Teuchos::RCP<const NOX::Abstract::MultiVector> X_imag;
  //     Teuchos::RCP<NOX::Abstract::MultiVector> U_real;
  //     Teuchos::RCP<NOX::Abstract::MultiVector> U_imag;
  //     const NOX::LAPACK::Vector *v1, *v2;
  //     NOX::LAPACK::Vector *w1, *w2;

  //     cX = dynamic_cast<const LOCA::Hopf::ComplexMultiVector*>(&X);
  //     X_real = cX->getRealMultiVec();
  //     X_imag = cX->getImagMultiVec();
  //     cU = dynamic_cast<LOCA::Hopf::ComplexMultiVector*>(&U);
  //     U_real = cU->getRealMultiVec();
  //     U_imag = cU->getImagMultiVec();

  //     // Concatenate X & Y into a single matrix
  //     NOX::LAPACK::Matrix< std::complex<double> > RHS(N,numCols);
  //     NOX::LAPACK::Matrix< std::complex<double> > LHS(N,numCols);
  //     for (int j=0; j<numCols; j++) {
  //       v1 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*X_real)[j]);
  //       v2 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*X_imag)[j]);
  //       for (int i=0; i<n; i++)
  //     RHS(i,j) = std::complex<double>((*v1)(i), (*v2)(i));
  //       for (int i=0; i<m; i++)
  //     RHS(i+n,j) = std::complex<double>(Y(i,j), Y(i+m,j));
  //     }

  //     // Solve for LHS
  //     augComplexSolver->apply(false, numCols, &RHS(0,0), &LHS(0,0));

  //     // Copy result into U and V
  //     for (int j=0; j<numCols; j++) {
  //       w1 = dynamic_cast<NOX::LAPACK::Vector*>(&(*U_real)[j]);
  //       w2 = dynamic_cast<NOX::LAPACK::Vector*>(&(*U_imag)[j]);
  //       for (int i=0; i<n; i++) {
  //     (*w1)(i) = LHS(i,j).real();
  //     (*w2)(i) = LHS(i,j).imag();
  //       }
  //       for (int i=0; i<m; i++) {
  //     V(i,j)   = LHS(n+i,j).real();
  //     V(i+m,j) = LHS(n+i,j).imag();
  //       }
  //     }
  //   }

  //   return NOX::Abstract::Group::Ok;

  // Compute J*X
  NOX::Abstract::Group::ReturnType status = op->apply(X, U);

  // Compute J*X + A*Y
  if (!isZeroA)
    U.update(Teuchos::NO_TRANS, 1.0, *A, Y, 1.0);

  // Compute B^T*X
  if (!isZeroB)
    B->multiplyDX(1.0, X, V);

  // Compute B^T*X + C*Y
  if (!isZeroC) {
    int e;
    if (isZeroB)
      e = V.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 0.0);
    else
      e = V.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 1.0);
    if (e < 0)
      status = NOX::Abstract::Group::Failed;
  }

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::applyTranspose(
    const NOX::Abstract::MultiVector& X,
    const NOX::Abstract::MultiVector::DenseMatrix& Y,
    NOX::Abstract::MultiVector& U,
    NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  // Compute J*X
  NOX::Abstract::Group::ReturnType status = op->applyTranspose(X, U);

  // Compute J*X + B*Y
  if (!isZeroA)
    B->addDX(Teuchos::NO_TRANS, 1.0, Y, 1.0, U);

  // Compute A^T*X
  if (!isZeroB)
    X.multiply(1.0, *A, V);

  // Compute A^T*X + C^T*Y
  if (!isZeroC) {
    int e;
    if (isZeroB)
      e = V.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 0.0);
    else
      e = V.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, *C, Y, 1.0);
    if (e < 0)
      status = NOX::Abstract::Group::Failed;
  }

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::applyInverse(
    Teuchos::ParameterList& params,
    const NOX::Abstract::MultiVector* F,
    const NOX::Abstract::MultiVector::DenseMatrix* G,
    NOX::Abstract::MultiVector& X,
    NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  return solve(false, params, F, G, X, Y);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::applyInverseTranspose(
    Teuchos::ParameterList& params,
    const NOX::Abstract::MultiVector* F,
    const NOX::Abstract::MultiVector::DenseMatrix* G,
    NOX::Abstract::MultiVector& X,
    NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  return solve(true, params, F, G, X, Y);
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::LAPACKDirectSolve::solve(
    bool trans,
    Teuchos::ParameterList& /* params */,
    const NOX::Abstract::MultiVector* F,
    const NOX::Abstract::MultiVector::DenseMatrix* G,
    NOX::Abstract::MultiVector& X,
    NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  bool isZeroF2 = (F == NULL);
  bool isZeroG2 = (G == NULL);

  // If F & G are zero, the solution is zero
  if (isZeroF2 && isZeroG2) {
    X.init(0.0);
    Y.putScalar(0.0);
    return NOX::Abstract::Group::Ok;
  }

  int numColsRHS;
  if (!isZeroF2)
    numColsRHS = F->numVectors();
  else
    numColsRHS = G->numCols();

  bool res = false;
  if (!isComplex) {
    const NOX::LAPACK::Vector *v;
    NOX::LAPACK::Vector *w;

    // Concatenate F & G into a single matrix
    NOX::LAPACK::Matrix<double> RHS(N,numColsRHS);
    if (isZeroF2) {
      for (int j=0; j<numColsRHS; j++)
        for (int i=0; i<n; i++)
          RHS(i,j) = 0.0;
    }
    else {
      for (int j=0; j<numColsRHS; j++) {
        v = dynamic_cast<const NOX::LAPACK::Vector*>(&(*F)[j]);
        TEUCHOS_ASSERT(v != NULL);
        for (int i=0; i<n; i++)
          RHS(i,j) = (*v)(i);
      }
    }
    if (isZeroG2) {
      for (int j=0; j<numColsRHS; j++)
        for (int i=0; i<m; i++)
          RHS(i+n,j) = 0.0;
    }
    else {
      for (int j=0; j<numColsRHS; j++)
        for (int i=0; i<m; i++)
          RHS(i+n,j) = (*G)(i,j);
    }

    // Solve for LHS
    res = augJacSolver->solve(trans, numColsRHS, &RHS(0,0));

    // Copy result into X and Y
    for (int j=0; j<numColsRHS; j++) {
      w = dynamic_cast<NOX::LAPACK::Vector*>(&X[j]);
      TEUCHOS_ASSERT(w != NULL);
      for (int i=0; i<n; i++)
        (*w)(i) = RHS(i,j);
      for (int i=0; i<m; i++)
        Y(i,j) = RHS(n+i,j);
    }
  }

  else {
#ifdef HAVE_TEUCHOS_COMPLEX
    const LOCA::Hopf::ComplexMultiVector* cF;
    LOCA::Hopf::ComplexMultiVector* cX;
    Teuchos::RCP<const NOX::Abstract::MultiVector> F_real;
    Teuchos::RCP<const NOX::Abstract::MultiVector> F_imag;
    Teuchos::RCP<NOX::Abstract::MultiVector> X_real;
    Teuchos::RCP<NOX::Abstract::MultiVector> X_imag;
    const NOX::LAPACK::Vector *v1, *v2;
    NOX::LAPACK::Vector *w1, *w2;

    if (!isZeroF2) {
      cF = dynamic_cast<const LOCA::Hopf::ComplexMultiVector*>(F);
      TEUCHOS_ASSERT(cF != NULL);
      F_real = cF->getRealMultiVec();
      F_imag = cF->getImagMultiVec();
    }
    cX = dynamic_cast<LOCA::Hopf::ComplexMultiVector*>(&X);
    TEUCHOS_ASSERT(cX != NULL);
    X_real = cX->getRealMultiVec();
    X_imag = cX->getImagMultiVec();

    // Concatenate F & G into a single matrix
    NOX::LAPACK::Matrix< std::complex<double> > RHS(N,numColsRHS);
    if (isZeroF2) {
      for (int j=0; j<numColsRHS; j++)
        for (int i=0; i<n; i++)
          RHS(i,j) = 0.0;
    }
    else {
      for (int j=0; j<numColsRHS; j++) {
        v1 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*F_real)[j]);
        TEUCHOS_ASSERT(v1 != NULL);
        v2 = dynamic_cast<const NOX::LAPACK::Vector*>(&(*F_imag)[j]);
        TEUCHOS_ASSERT(v2 != NULL);
        for (int i=0; i<n; i++)
          RHS(i,j) = std::complex<double>((*v1)(i), (*v2)(i));
      }
    }
    if (isZeroG2) {
      for (int j=0; j<numColsRHS; j++)
        for (int i=0; i<m; i++)
          RHS(i+n,j) = 0.0;
    }
    else {
      for (int j=0; j<numColsRHS; j++)
        for (int i=0; i<m; i++)
          RHS(i+n,j) = std::complex<double>((*G)(i,j), (*G)(i+m,j));
    }

    // Solve for LHS
    res = augComplexSolver->solve(trans, numColsRHS, &RHS(0,0));

    // Copy result into X and Y
    for (int j=0; j<numColsRHS; j++) {
      w1 = dynamic_cast<NOX::LAPACK::Vector*>(&(*X_real)[j]);
      TEUCHOS_ASSERT(w1 != NULL);
      w2 = dynamic_cast<NOX::LAPACK::Vector*>(&(*X_imag)[j]);
      TEUCHOS_ASSERT(w2 != NULL);
      for (int i=0; i<n; i++) {
        (*w1)(i) = RHS(i,j).real();
        (*w2)(i) = RHS(i,j).imag();
      }
      for (int i=0; i<m; i++) {
        Y(i,j)   = RHS(n+i,j).real();
        Y(i+m,j) = RHS(n+i,j).imag();
      }
    }
#else
    globalData->locaErrorCheck->throwError(
        "LOCA::BorderedSolver::LAPACKDirectSolve::solve()",
        "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
#endif
  }

  if (!res)
    return NOX::Abstract::Group::Failed;
  else
    return NOX::Abstract::Group::Ok;
}
