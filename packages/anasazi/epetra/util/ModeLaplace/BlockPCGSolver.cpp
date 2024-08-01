// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#include "BlockPCGSolver.h"


BlockPCGSolver::BlockPCGSolver(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                               double _tol, int _iMax, int _verb)
               : MyComm(_Comm),
                 callBLAS(),
                 callLAPACK(),
                 K(KK),
                 Prec(0),
                 tolCG(_tol),
                 iterMax(_iMax),
                 verbose(_verb),
                 workSpace(0),
                 lWorkSpace(0),
                 numSolve(0),
                 maxIter(0),
                 sumIter(0),
                 minIter(10000)
               {

}


BlockPCGSolver::BlockPCGSolver(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                               Epetra_Operator *PP,
                               double _tol, int _iMax, int _verb)
               : MyComm(_Comm),
                 callBLAS(),
                 callLAPACK(),
                 K(KK),
                 Prec(PP),
                 tolCG(_tol),
                 iterMax(_iMax),
                 verbose(_verb),
                 workSpace(0),
                 lWorkSpace(0),
                 numSolve(0),
                 maxIter(0),
                 sumIter(0),
                 minIter(10000)
               {}


BlockPCGSolver::~BlockPCGSolver() {
  if (workSpace) {
    delete[] workSpace;
    workSpace = 0;
    lWorkSpace = 0;
  }
}


void BlockPCGSolver::setPreconditioner(Epetra_Operator *PP) {
  Prec = PP;
}


int BlockPCGSolver::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  return K->Apply(X, Y);
}


int BlockPCGSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  int xcol = X.NumVectors();
  int info = 0;

  if (Y.NumVectors() < xcol)
    return -1;

  // Use block PCG for multiple right-hand sides
  info = (xcol == 1) ? Solve(X, Y) : Solve(X, Y, xcol);

  return info;
}


int BlockPCGSolver::Solve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const {

  int info = 0;
  int localVerbose = verbose*(MyComm.MyPID() == 0);

  int xr = X.MyLength();

  int wSize = 3*xr;

  if (lWorkSpace < wSize) {
    if (workSpace)
      delete[] workSpace;
    workSpace = new (std::nothrow) double[wSize];
    if (workSpace == 0) {
      info = -1;
      return info;
    }
    lWorkSpace = wSize;
  } // if (lWorkSpace < wSize)

  double *pointer = workSpace;

  Epetra_Vector r(Epetra_DataAccess::View, X.Map(), pointer);
  pointer = pointer + xr;

  Epetra_Vector p(Epetra_DataAccess::View, X.Map(), pointer);
  pointer = pointer + xr;

  // Note: Kp and z uses the same memory space
  Epetra_Vector Kp(Epetra_DataAccess::View, X.Map(), pointer);
  Epetra_Vector z(Epetra_DataAccess::View, X.Map(), pointer);

  double tmp;
  double initNorm = 0.0, rNorm = 0.0, newRZ = 0.0, oldRZ = 0.0, alpha = 0.0;
  double tolSquare = tolCG*tolCG;

  memcpy(r.Values(), X.Values(), xr*sizeof(double));
  tmp = callBLAS.DOT(xr, r.Values(), 1, r.Values(), 1);
  MyComm.SumAll(&tmp, &initNorm, 1);

  Y.PutScalar(0.0);

  if (localVerbose > 1) {
    std::cout << std::endl;
    std::cout  << " --- PCG Iterations --- " << std::endl;
  }

  int iter;
  for (iter = 1; iter <= iterMax; ++iter) {

    if (Prec) {
      Prec->ApplyInverse(r, z);
    }
    else {
      memcpy(z.Values(), r.Values(), xr*sizeof(double));
    }

    if (iter == 1) {
      tmp = callBLAS.DOT(xr, r.Values(), 1, z.Values(), 1);
      MyComm.SumAll(&tmp, &newRZ, 1);
      memcpy(p.Values(), z.Values(), xr*sizeof(double));
    }
    else {
      oldRZ = newRZ;
      tmp = callBLAS.DOT(xr, r.Values(), 1, z.Values(), 1);
      MyComm.SumAll(&tmp, &newRZ, 1);
      p.Update(1.0, z, newRZ/oldRZ);
    }

    K->Apply(p, Kp);

    tmp = callBLAS.DOT(xr, p.Values(), 1, Kp.Values(), 1);
    MyComm.SumAll(&tmp, &alpha, 1);
    alpha = newRZ/alpha;

    TEUCHOS_TEST_FOR_EXCEPTION(alpha <= 0.0, std::runtime_error,
                         " !!! Non-positive value for p^TKp (" << alpha << ") !!!");

    callBLAS.AXPY(xr, alpha, p.Values(), 1, Y.Values(), 1);

    alpha *= -1.0;
    callBLAS.AXPY(xr, alpha, Kp.Values(), 1, r.Values(), 1);

    // Check convergence
    tmp = callBLAS.DOT(xr, r.Values(), 1, r.Values(), 1);
    MyComm.SumAll(&tmp, &rNorm, 1);

    if (localVerbose > 1) {
      std::cout  << "   Iter. " << iter;
      std::cout.precision(4);
      std::cout.setf(std::ios::scientific, std::ios::floatfield);
      std::cout << " Residual reduction " << std::sqrt(rNorm/initNorm) << std::endl;
    }

    if (rNorm <= tolSquare*initNorm)
      break;

  } // for (iter = 1; iter <= iterMax; ++iter)

  if (localVerbose == 1) {
    std::cout << std::endl;
    std::cout << " --- End of PCG solve ---" << std::endl;
    std::cout << "   Iter. " << iter;
    std::cout.precision(4);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << " Residual reduction " << std::sqrt(rNorm/initNorm) << std::endl;
    std::cout << std::endl;
  }

  if (localVerbose > 1) {
    std::cout << std::endl;
  }

  numSolve += 1;

  minIter = (iter < minIter) ? iter : minIter;
  maxIter = (iter > maxIter) ? iter : maxIter;
  sumIter += iter;

  return info;
}


int BlockPCGSolver::Solve(const Epetra_MultiVector &X, Epetra_MultiVector &Y, int blkSize) const {

  int xrow = X.MyLength();
  int xcol = X.NumVectors();
  int ycol = Y.NumVectors();

  int info = 0;
  int localVerbose = verbose*(MyComm.MyPID() == 0);
  double *valX = X.Values();
  int NB = 3 + callLAPACK.ILAENV(1, "hetrd", "u", blkSize);
  int lworkD = (blkSize > NB) ? blkSize*blkSize : NB*blkSize;
  int wSize = 4*blkSize*xrow + 3*blkSize + 2*blkSize*blkSize + lworkD;

  bool useY = true;
  if (ycol % blkSize != 0) {
    // Allocate an extra block to store the solutions
    wSize += blkSize*xrow;
    useY = false;
  }

  if (lWorkSpace < wSize) {
    delete[] workSpace;
    workSpace = new (std::nothrow) double[wSize];
    if (workSpace == 0) {
      info = -1;
      return info;
    }
    lWorkSpace = wSize;
  } // if (lWorkSpace < wSize)

  double *pointer = workSpace;

  // Array to store the matrix PtKP
  double *PtKP = pointer;
  pointer = pointer + blkSize*blkSize;

  // Array to store coefficient matrices
  double *coeff = pointer;
  pointer = pointer + blkSize*blkSize;

  // Workspace array
  double *workD = pointer;
  pointer = pointer + lworkD;

  // Array to store the eigenvalues of P^t K P
  double *da = pointer;
  pointer = pointer + blkSize;

  // Array to store the norms of right hand sides
  double *initNorm = pointer;
  pointer = pointer + blkSize;

  // Array to store the norms of residuals
  double *resNorm = pointer;
  pointer = pointer + blkSize;

  // Array to store the residuals
  double *valR = pointer;
  pointer = pointer + xrow*blkSize;
  Epetra_MultiVector R(Epetra_DataAccess::View, X.Map(), valR, xrow, blkSize);

  // Array to store the preconditioned residuals
  double *valZ = pointer;
  pointer = pointer + xrow*blkSize;
  Epetra_MultiVector Z(Epetra_DataAccess::View, X.Map(), valZ, xrow, blkSize);

  // Array to store the search directions
  double *valP = pointer;
  pointer = pointer + xrow*blkSize;
  Epetra_MultiVector P(Epetra_DataAccess::View, X.Map(), valP, xrow, blkSize);

  // Array to store the image of the search directions
  double *valKP = pointer;
  pointer = pointer + xrow*blkSize;
  Epetra_MultiVector KP(Epetra_DataAccess::View, X.Map(), valKP, xrow, blkSize);

  // Pointer to store the solutions
  double *valSOL = (useY == true) ? Y.Values() : pointer;

  int iRHS;
  for (iRHS = 0; iRHS < xcol; iRHS += blkSize) {

    int numVec = (iRHS + blkSize < xcol) ? blkSize : xcol - iRHS;

    // Set the initial residuals to the right hand sides
    if (numVec < blkSize) {
      R.Random();
    }
    memcpy(valR, valX + iRHS*xrow, numVec*xrow*sizeof(double));

    // Set the initial guess to zero
    valSOL = (useY == true) ? Y.Values() + iRHS*xrow : valSOL;
    Epetra_MultiVector SOL(Epetra_DataAccess::View, X.Map(), valSOL, xrow, blkSize);
    SOL.PutScalar(0.0);

    int ii = 0;
    int iter = 0;
    int nFound = 0;

    R.Norm2(initNorm);

    if (localVerbose > 1) {
      std::cout << std::endl;
      std::cout << " Vectors " << iRHS << " to " << iRHS + numVec - 1 << std::endl;
      if (localVerbose > 2) {
        std::fprintf(stderr,"\n");
        for (ii = 0; ii < numVec; ++ii) {
          std::cout << " ... Initial Residual Norm " << ii << " = " << initNorm[ii] << std::endl;
        }
        std::cout << std::endl;
      }
    }

    // Iteration loop
    for (iter = 1; iter <= iterMax; ++iter) {

      // Apply the preconditioner
      if (Prec)
        Prec->ApplyInverse(R, Z);
      else
        Z = R;

      // Define the new search directions
      if (iter == 1) {
        P = Z;
      }
      else {
        // Compute P^t K Z
        callBLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, blkSize, blkSize, xrow, 1.0, KP.Values(), xrow, Z.Values(), xrow,
                      0.0, workD, blkSize);
        MyComm.SumAll(workD, coeff, blkSize*blkSize);

        // Compute the coefficient (P^t K P)^{-1} P^t K Z
        callBLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, blkSize, blkSize, blkSize, 1.0, PtKP, blkSize, coeff, blkSize,
                      0.0, workD, blkSize);
        for (ii = 0; ii < blkSize; ++ii)
          callBLAS.SCAL(blkSize, da[ii], workD + ii, blkSize);
        callBLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, blkSize, blkSize, blkSize, 1.0, PtKP, blkSize, workD, blkSize,
                      0.0, coeff, blkSize);

        // Update the search directions
        // Note: Use KP as a workspace
        memcpy(KP.Values(), P.Values(), xrow*blkSize*sizeof(double));
        callBLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, xrow, blkSize, blkSize, 1.0, KP.Values(), xrow, coeff, blkSize,
                      0.0, P.Values(), xrow);

        P.Update(1.0, Z, -1.0);

      } // if (iter == 1)

      K->Apply(P, KP);

      // Compute P^t K P
      callBLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, blkSize, blkSize, xrow, 1.0, P.Values(), xrow, KP.Values(), xrow,
                    0.0, workD, blkSize);
      MyComm.SumAll(workD, PtKP, blkSize*blkSize);

      // Eigenvalue decomposition of P^t K P
      callLAPACK.SYEV('V', 'U', blkSize, PtKP, blkSize, da, workD, lworkD, &info);
      if (info) {
        // Break the loop as spectral decomposition failed
        break;
      } // if (info)

      // Compute the pseudo-inverse of the eigenvalues
      for (ii = 0; ii < blkSize; ++ii) {
        TEUCHOS_TEST_FOR_EXCEPTION(da[ii] < 0.0, std::runtime_error, "Negative "
                           "eigenvalue for P^T K P: da[" << ii << "] = "
                           << da[ii] << ".");
        da[ii] = (da[ii] == 0.0) ? 0.0 : 1.0/da[ii];
      } // for (ii = 0; ii < blkSize; ++ii)

      // Compute P^t R
      callBLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, blkSize, blkSize, xrow, 1.0, P.Values(), xrow, R.Values(), xrow,
                    0.0, workD, blkSize);
      MyComm.SumAll(workD, coeff, blkSize*blkSize);

      // Compute the coefficient (P^t K P)^{-1} P^t R
      callBLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, blkSize, blkSize, blkSize, 1.0, PtKP, blkSize, coeff, blkSize,
                    0.0, workD, blkSize);
      for (ii = 0; ii < blkSize; ++ii)
        callBLAS.SCAL(blkSize, da[ii], workD + ii, blkSize);
      callBLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, blkSize, blkSize, blkSize, 1.0, PtKP, blkSize, workD, blkSize,
                    0.0, coeff, blkSize);

      // Update the solutions
      callBLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, xrow, blkSize, blkSize, 1.0, P.Values(), xrow, coeff, blkSize,
                    1.0, valSOL, xrow);

      // Update the residuals
      callBLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, xrow, blkSize, blkSize, -1.0, KP.Values(), xrow, coeff, blkSize,
                    1.0, R.Values(), xrow);

      // Check convergence
      R.Norm2(resNorm);
      nFound = 0;
      for (ii = 0; ii < numVec; ++ii) {
        if (resNorm[ii] <= tolCG*initNorm[ii])
          nFound += 1;
      }

      if (localVerbose > 1) {
        std::cout << " Vectors " << iRHS << " to " << iRHS + numVec - 1;
        std::cout << " -- Iteration " << iter << " -- " << nFound << " converged vectors\n";
        if (localVerbose > 2) {
          std::cout << std::endl;
          for (ii = 0; ii < numVec; ++ii) {
            std::cout << " ... ";
            std::cout.width(5);
            std::cout << ii << " ... Residual = ";
            std::cout.precision(2);
            std::cout.setf(std::ios::scientific, std::ios::floatfield);
            std::cout << resNorm[ii] << " ... Right Hand Side = " << initNorm[ii] << std::endl;
          }
          std::cout << std::endl;
        }
      }

      if (nFound == numVec) {
        break;
      }

    }  // for (iter = 1; iter <= maxIter; ++iter)

    if (useY == false) {
      // Copy the solutions back into Y
      memcpy(Y.Values() + xrow*iRHS, valSOL, numVec*xrow*sizeof(double));
    }

    numSolve += nFound;

    if (nFound == numVec) {
      minIter = (iter < minIter) ? iter : minIter;
      maxIter = (iter > maxIter) ? iter : maxIter;
      sumIter += iter;
    }

  } // for (iRHS = 0; iRHS < xcol; iRHS += blkSize)

  return info;
}


