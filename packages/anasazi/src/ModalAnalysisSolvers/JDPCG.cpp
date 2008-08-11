//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the cout
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical cout SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://software.sandia.gov/trilinos/ ).
//
// The distribution of this software follows also the rules defined in Trilinos.
// This notice shall be marked on any reproduction of this software, in whole or
// in part.
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#include "JDPCG.h"


JDPCG::JDPCG(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     const Epetra_Operator *MM, const Epetra_Operator *PP, int _blk, int _numBlk,
                     double _tol, int _maxIterES, int _maxIterLS, int _verb, double *_weight) :
           myVerify(_Comm),
           callBLAS(),
           callFortran(),
           callLAPACK(),
           modalTool(_Comm),
           mySort(),
           MyComm(_Comm),
           K(KK),
           M(MM),
           Prec(PP),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIterES),
           maxIterLinearSolve(_maxIterLS),
           blockSize(_blk),
           numBlock(_numBlk),
           normWeight(_weight),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           maxSpaceSize(0),
           sumSpaceSize(0),
           spaceSizeHistory(0),
           maxIterPCG(0),
           sumIterPCG(0),
           iterPCGHistory(0),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           numCorrectionPrec(0),
           numCorrectionSolve(0),
           numPCGmassOp(0),
           numPCGstifOp(0),
           numRestart(0),
           outerIter(0),
           precOp(0),
           residual(0),
           stifOp(0),
           timeBuildQtMPMQ(0.0),
           timeCorrectionPrec(0.0),
           timeCorrectionSolve(0.0),
           timeLocalProj(0.0),
           timeLocalSolve(0.0),
           timeLocalUpdate(0.0),
           timeMassOp(0.0),
           timeNorm(0.0),
           timeOrtho(0.0),
           timeOuterLoop(0.0),
           timePCGEigCheck(0.0),
           timePCGLoop(0.0),
           timePCGOpMult(0.0),
           timePCGPrec(0.0),
           timePostProce(0.0),
           timePrecOp(0.0),
           timeResidual(0.0),
           timeRestart(0.0),
           timeStifOp(0.0)
           {

}


JDPCG::~JDPCG() {

  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

  if (spaceSizeHistory) {
    delete[] spaceSizeHistory;
    spaceSizeHistory = 0;
  }

  if (iterPCGHistory) {
    delete[] iterPCGHistory;
    iterPCGHistory = 0;
  }

}


int JDPCG::jacobiPreconditioner(const Epetra_MultiVector &B, Epetra_MultiVector &PrecB,
                       const Epetra_MultiVector *U, const Epetra_MultiVector *V,
                       double *invQtMPMQ, int ldQtMPMQ, double *PMQ, double *work, double *WS) {

  // This routine applies the "projected" preconditioner to the vectors B and stores
  // the result in the vectors PrecB.
  //
  //             PrecB = (I - P^{-1}MQ ( Q^tMP^{-1}MQ )^{-1} Q^t M) P^{-1} * B
  //
  // where P is the preconditioner.
  //
  // B = Input vectors
  // PrecB = Preconditioned vectors
  // V = converged eigenvectors
  // U = tentative eigenvectors to be corrected.
  //
  // Note: Q = [V, U]
  //
  // invQtMPMQ = Factor form of the matrix Q^t*M*P^{-1}*M*Q (from POTRF)
  // ldQtMPMQ = Leading dimension for invQtMPMQ
  //
  // PMQ = Array to store the column vectors P^{-1}*M*Q
  //      Assumption: PMQ has the same distribution than B over the processors
  //
  // work = Workspace array of size 2 * (# of columns in Q) * (# of columns in B) 
  //
  // WS = Workspace array of size (# of rows in B) * (# of columns in B)

  int info = 0;

  int bC = B.NumVectors();
  int bR = B.MyLength();

  // Compute P^{-1} B
  timeCorrectionPrec -= MyWatch.WallTime();
  if (Prec) {
    Prec->ApplyInverse(B, PrecB);
  }
  else {
    memcpy(PrecB.Values(), B.Values(), bR*bC*sizeof(double));
  }
  timeCorrectionPrec += MyWatch.WallTime();

  int uC = U->NumVectors();
  int vC = (V) ? V->NumVectors() : 0;

  int sizeQ = uC + vC;

  // Compute Q^t M P^{-t} B
  callBLAS.GEMM('T', 'N', sizeQ, bC, bR, 1.0, PMQ, bR, B.Values(), bR, 
                0.0, work + sizeQ*bC, sizeQ);
  MyComm.SumAll(work + sizeQ*bC, work, sizeQ*bC);

  memcpy(work + sizeQ*bC, work, sizeQ*bC*sizeof(double));

  // Compute ( Q^t M P^{-1} M Q )^{-1} Q^t M P^{-t} B
  callLAPACK.POTRS('U', sizeQ, bC, invQtMPMQ, ldQtMPMQ, work, sizeQ, &info);

  // Subtract  P^{-1} M Q ( Q^t M P^{-1} M Q )^{-1} Q^t M P^{-t} B from P^{-t} B
  callBLAS.GEMM('N', 'N', bR, bC, sizeQ, -1.0, PMQ, bR, work, sizeQ, 1.0, PrecB.Values(), bR);

  Epetra_MultiVector MPB(View, B.Map(), WS, bR, bC);
  M->Apply(PrecB, MPB);

  // Compute Q^t M PrecB
  callBLAS.GEMM('T', 'N', uC, bC, bR, 1.0, U->Values(), bR, MPB.Values(), bR, 
                0.0, work + vC + sizeQ*bC, sizeQ);
  if (V) {
    callBLAS.GEMM('T', 'N', vC, bC, bR, 1.0, V->Values(), bR, MPB.Values(), bR, 
                  0.0, work + sizeQ*bC, sizeQ);
  }
  MyComm.SumAll(work + sizeQ*bC, work, sizeQ*bC);

  // Test for second projection
  double *Mnorm = new double[bC];
  MPB.Dot(PrecB, Mnorm);

  bool doSecond = false;
  for (int j = 0; j < bC; ++j) {
    double tolOrtho = 1.0e-28*Mnorm[j];
    for (int i = 0; i < sizeQ; ++i) {
      double tmp = work[i + j*sizeQ];
      if (tmp*tmp > tolOrtho) {
        doSecond = true;
        break;
      }
    } 
    if (doSecond == true) {
      // Compute ( Q^t M P^{-1} M Q )^{-1} Q^t M PrecB
      callLAPACK.POTRS('U', sizeQ, bC, invQtMPMQ, ldQtMPMQ, work, sizeQ, &info);
      // Subtract  P^{-1} M Q ( Q^t M P^{-1} M Q )^{-1} Q^t M PrecB from PrecB
      callBLAS.GEMM('N', 'N', bR, bC, sizeQ, -1.0, PMQ, bR, work, sizeQ, 1.0, PrecB.Values(), bR);
      break;
    }
  }
  delete[] Mnorm;

  numCorrectionPrec += bC;

  return info;

}


int JDPCG::jacobiPCG(Epetra_MultiVector &Rlin, Epetra_MultiVector &Y,
                     const Epetra_MultiVector *U, const Epetra_MultiVector *V, 
                     double eta, double tolCG, int iterMax, 
                     double *invQtMPMQ, int ldQtMPMQ, double *PMQ,
                     double *workPrec, double *workPCG,
                     const Epetra_Vector *vectWeight, const Epetra_MultiVector *orthoVec) {

  // This routine applies a block PCG algorithm to solve the equation
  //
  //      (I - MQ*Q^t) * ( K - eta * M ) * (I - Q*Q^t*M) Y = X
  //
  // with (I - Q*Q^t*M) * Y = Y
  // with Q = [V, U]
  // where the preconditioner is given by
  //
  //      (I - MQ*Q^t) * Prec^{-1} * (I - Q*Q^t*M)
  //
  // Rlin = Input vectors
  // Y = Solution vectors
  // V = converged eigenvectors
  // U = tentative eigenvectors to be corrected.
  //
  // eta = shift for the linear operator
  //
  // tolCG = Tolerance required for convergence
  // iterMax = Maximum number of iterations allowed
  //
  // invQtMPMQ = Factor form of the matrix Q^t*M*P^{-1}*M*Q (from POTRF)
  // ldQtMPMQ = Leading dimension for invQtMPMQ
  //
  // PMQ = Array to store the column vectors P^{-1}*M*Q
  //      Assumption: PMQ has the same distribution than B over the processors
  //
  // workPrec = Workspace array of size 2 * (# of columns in Q) * (# of columns in X)
  //            This workspace is exclusively used to apply the preconditioner
  //
  // workPCG = Workspace array for the variables in the PCG algorithm
  //           Its size must allow the definition of
  //           - 5 blocks of (# of columns in X) vectors distributed as X across the processors
  //           - 4 arrays of length (# of columns in X)
  //           - 3 square matrices of size (# of columns in X)
  //
  // vectWeight = Weights for the L^2 norm to compute to check the eigenresiduals
  //
  // orthoVec = Space where CG computations are orthogonal to.

  int xrow = Y.MyLength();
  int xcol = Y.NumVectors();

  int info = 0;
  int localVerbose = verbose*(MyComm.MyPID() == 0);

  double *pointer = workPCG;

  // Arrays associated with the solution to the linear system

  // Array to store the matrix PtKP
  double *PtKP = pointer;
  pointer = pointer + xcol*xcol;

  // Array to store coefficient matrices
  double *coeff = pointer;
  pointer = pointer + xcol*xcol;

  // Workspace array
  double *workD = pointer;
  pointer = pointer + xcol*xcol;

  // Array to store the eigenvalues of P^t K P
  double *da = pointer;
  pointer = pointer + xcol;

  // Array to store the norms of right hand sides
  double *initNorm = pointer;
  pointer = pointer + xcol;

  // Array to store the norms of current residuals
  double *resNorm = pointer;
  pointer = pointer + xcol;

  // Array to store the preconditioned residuals
  double *valZ = pointer;
  pointer = pointer + xrow*xcol;
  Epetra_MultiVector Z(View, Y.Map(), valZ, xrow, xcol);

  // Array to store the search directions
  double *valP = pointer;
  pointer = pointer + xrow*xcol;
  Epetra_MultiVector P(View, Y.Map(), valP, xrow, xcol);

  // Array to store the image of the search directions
  double *valKP = pointer;
  pointer = pointer + xrow*xcol;
  Epetra_MultiVector KP(View, Y.Map(), valKP, xrow, xcol);

  // Arrays associated to the corrected eigenvectors
  // Array to store the projected stiffness matrix
  double *UtKU = pointer;
  pointer = pointer + xcol*xcol;
  
  // Array to store the corrected eigenvalues
  double *theta = pointer;
  pointer = pointer + xcol;
  
  // Array to store the norms of eigen-residuals for corrected vectors
  double *resEig = pointer;
  pointer = pointer + xcol;
  
  // Array to store the norms of previous eigen-residuals for corrected vectors
  double *oldEig = pointer;
  pointer = pointer + xcol;
  
  // Array to store the stiffness-image of the corrected eigenvectors
  double *valKU = pointer;
  pointer = pointer + xrow*xcol;
  Epetra_MultiVector KU(View, Y.Map(), valKU, xrow, xcol);

  // Array to store the mass-image of the corrected eigenvectors
  Epetra_MultiVector MU(View, Y.Map(), (M) ? pointer : Z.Values(), xrow, xcol);

  // Set the initial guess to zero
  Y.PutScalar(0.0);

  int ii;
  int iter;
  int nFound;

  // Define the size of the "orthogonal" space [V, U]
  // Note: We assume U is non zero and thus sizeQ is non zero.
  int sizeQ = 0;
  sizeQ += (U) ? U->NumVectors() : 0;
  sizeQ += (V) ? V->NumVectors() : 0;

  bool isNegative = false;

  Rlin.Norm2(initNorm);

  if (localVerbose > 3) {
    cout << endl;
    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);
    for (ii = 0; ii < xcol; ++ii) {
      cout << " ... Initial Residual Norm " << ii << " = " << initNorm[ii] << endl;
    }
    cout << endl;
  }

  // Iteration loop
  timePCGLoop -= MyWatch.WallTime();
  for (iter = 1; iter <= iterMax; ++iter) {

    // Apply the preconditioner
    timePCGPrec -= MyWatch.WallTime();
    jacobiPreconditioner(Rlin, Z, U, V, invQtMPMQ, ldQtMPMQ, PMQ, workPrec, MU.Values());
    timePCGPrec += MyWatch.WallTime();

    if (orthoVec) {
      // Note: Use MU as workspace
      if (M)
        M->Apply(Z, MU);
      modalTool.massOrthonormalize(Z, MU, M, *orthoVec, blockSize, 1);
    }

    // Define the new search directions
    if (iter == 1) {
      P = Z;
    }
    else {
      // Compute P^t K Z
      callBLAS.GEMM('T', 'N', xcol, xcol, xrow, 1.0, KP.Values(), xrow, Z.Values(), xrow,
                    0.0, workD, xcol);
      MyComm.SumAll(workD, coeff, xcol*xcol);

      // Compute the coefficient (P^t K P)^{-1} P^t K Z
      callBLAS.GEMM('T', 'N', xcol, xcol, xcol, 1.0, PtKP, xcol, coeff, xcol,
                    0.0, workD, xcol);
      for (ii = 0; ii < xcol; ++ii)
        callFortran.SCAL_INCX(xcol, da[ii], workD + ii, xcol);
      callBLAS.GEMM('N', 'N', xcol, xcol, xcol, 1.0, PtKP, xcol, workD, xcol,
                    0.0, coeff, xcol);

      // Update the search directions 
      // Note: Use KP as a workspace
      memcpy(KP.Values(), P.Values(), xrow*xcol*sizeof(double));
      callBLAS.GEMM('N', 'N', xrow, xcol, xcol, 1.0, KP.Values(), xrow, coeff, xcol,
                    0.0, P.Values(), xrow);

      P.Update(1.0, Z, -1.0);

    } // if (iter == 1)

    timePCGOpMult -= MyWatch.WallTime();
    K->Apply(P, KP);
    numPCGstifOp += xcol;
    if (eta != 0.0) {
      // Apply the mass matrix
      // Note: Use Z as a workspace
      if (M) {
        M->Apply(P, Z);
        callBLAS.AXPY(xrow*xcol, -eta, Z.Values(), KP.Values());
      }
      else {
        callBLAS.AXPY(xrow*xcol, -eta, P.Values(), KP.Values());
      }
      numPCGmassOp += xcol;
    }
    timePCGOpMult += MyWatch.WallTime();

    // Compute P^t K P
    callBLAS.GEMM('T', 'N', xcol, xcol, xrow, 1.0, P.Values(), xrow, KP.Values(), xrow, 
                  0.0, workD, xcol);
    MyComm.SumAll(workD, PtKP, xcol*xcol);

    // Eigenvalue decomposition of P^t K P
    int nev = xcol;
    info = modalTool.directSolver(xcol, PtKP, xcol, 0, 0, nev, PtKP, xcol, da, 0, 10);

    if (info != 0) {
      // Break the loop as spectral decomposition failed
      info = - iterMax - 1;
      sumIterPCG += iter;
      maxIterPCG = (iter > maxIterPCG) ? iter : maxIterPCG;
      return info;
    } // if (info != 0)

    // Compute the pseudo-inverse of the eigenvalues
    for (ii = 0; ii < xcol; ++ii) {
      if (da[ii] < 0.0) {
        isNegative = true;
        break;
      }
      else {
        da[ii] = (da[ii] == 0.0) ? 0.0 : 1.0/da[ii];
      }
    } // for (ii = 0; ii < xcol; ++ii)

    if (isNegative == true) {
      if (localVerbose > 0) {
        cout << endl;
        cout.precision(4);
        cout.setf(ios::scientific, ios::floatfield);
        cout << " !! Negative eigenvalue in block PCG (" << da[ii] << ") !!\n";
        cout << endl;
      }
      info = - iter;
      sumIterPCG += iter;
      maxIterPCG = (iter > maxIterPCG) ? iter : maxIterPCG;
      return info;
    }

    // Compute P^t R
    callBLAS.GEMM('T', 'N', xcol, xcol, xrow, 1.0, P.Values(), xrow, Rlin.Values(), xrow,
                  0.0, workD, xcol);
    MyComm.SumAll(workD, coeff, xcol*xcol);

    // Compute the coefficient (P^t K P)^{-1} P^t R
    callBLAS.GEMM('T', 'N', xcol, xcol, xcol, 1.0, PtKP, xcol, coeff, xcol,
                  0.0, workD, xcol);
    for (ii = 0; ii < xcol; ++ii)
      callFortran.SCAL_INCX(xcol, da[ii], workD + ii, xcol);
    callBLAS.GEMM('N', 'N', xcol, xcol, xcol, 1.0, PtKP, xcol, workD, xcol,
                  0.0, coeff, xcol);

    // Update the solutions of the linear system
    callBLAS.GEMM('N', 'N', xrow, xcol, xcol, 1.0, P.Values(), xrow, coeff, xcol,
                  1.0, Y.Values(), xrow);

    // Update the residuals for the linear system
    callBLAS.GEMM('N', 'N', xrow, xcol, xcol, -1.0, KP.Values(), xrow, coeff, xcol,
                  1.0, Rlin.Values(), xrow);

    // Check convergence 
    Rlin.Norm2(resNorm);
    nFound = 0;
    for (ii = 0; ii < xcol; ++ii) {
      if (resNorm[ii] <= tolCG*initNorm[ii])
        nFound += 1;
    }

    if (localVerbose > 3) {
      cout << endl;
      for (ii = 0; ii < xcol; ++ii) {
        cout << " ... ";
        cout.width(5);
        cout << ii << " ... Residual = ";
        cout.precision(4);
        cout.setf(ios::scientific, ios::floatfield);
        cout << resNorm[ii] << " ... Right Hand Side = " << initNorm[ii] << endl;
      }
      cout << endl;
    }

    if (nFound == xcol) {
      info = iter;
      break;
    }

    timePCGEigCheck -= MyWatch.WallTime();
    
    // Check the residuals for the corrected eigenvectors
    // Note: Use Z as workspace to store the corrected vectors
    Z.Update(1.0, *U, -1.0, Y, 0.0);

    // Compute U^t K U
    K->Apply(Z, KU);
    numPCGstifOp += xcol;
    callBLAS.GEMM('T', 'N', xcol, xcol, xrow, 1.0, Z.Values(), xrow, KU.Values(), xrow,
                  0.0, workD, xcol);
    MyComm.SumAll(workD, UtKU, xcol*xcol);

    // Compute U^t M U
    // Note: Use coeff as storage space
    if (M) {
      M->Apply(Z, MU);
      numPCGmassOp += xcol;
      callBLAS.GEMM('T', 'N', xcol, xcol, xrow, 1.0, Z.Values(), xrow, MU.Values(), xrow,
                    0.0, workD, xcol);
      MyComm.SumAll(workD, coeff, xcol*xcol);
    }
    else {
      callBLAS.GEMM('T', 'N', xcol, xcol, xrow, 1.0, Z.Values(), xrow, Z.Values(), xrow,
                    0.0, workD, xcol);
      MyComm.SumAll(workD, coeff, xcol*xcol);
    }

    nev = xcol;
    info = modalTool.directSolver(xcol, UtKU, xcol, coeff, xcol, nev, workD, xcol, theta, 0, 
                                  (blockSize == 1) ? 1 : 0);

    if ((info < 0) || (nev < xcol)) {
      // Break the loop as spectral decomposition failed
      info = - iterMax - 1;
      sumIterPCG += iter;
      maxIterPCG = (iter > maxIterPCG) ? iter : maxIterPCG;
      return info; 
    } // if ((info < 0) || (nev < xcol))
    
    // Compute the eigenresiduals for the corrected vectors
    // Note: Use Z as workspace to store the residuals
    if (M) {
      callBLAS.GEMM('N', 'N', xrow, xcol, xcol, 1.0, MU.Values(), xrow, workD, xcol,
                    0.0, Z.Values(), xrow);
    }
    for (ii = 0; ii < xcol; ++ii)
      callBLAS.SCAL(xrow, theta[ii], Z.Values() + ii*xrow);
    callBLAS.GEMM('N', 'N', xrow, xcol, xcol, 1.0, KU.Values(), xrow, workD, xcol,
                  -1.0, Z.Values(), xrow);
    
    if (vectWeight)
      Z.NormWeighted(*vectWeight, resEig);
    else
      Z.Norm2(resEig);
    
    timePCGEigCheck += MyWatch.WallTime();
    
    if (iter > 1) {
      // Scale the norms of residuals with the eigenvalues
      // Count the number of converged eigenvectors
      nFound = 0;
      int nGrow = 0;
      for (ii = 0; ii < xcol; ++ii) {
        nFound = (resEig[ii] < tolEigenSolve*theta[ii]) ? nFound + 1 : nFound;
        nGrow = (resEig[ii] > oldEig[ii]) ? nGrow + 1 : nGrow;
      } // for (ii = 0; ii < xcol; ++ii)
      if ((nFound > 0) || (nGrow > 0)) {
        info = iter;
        break;
      }
    } // if (iter > 1)
    
    memcpy(oldEig, resEig, xcol*sizeof(double));
    
  }  // for (iter = 1; iter <= maxIter; ++iter)
  timePCGLoop += MyWatch.WallTime();

  sumIterPCG += iter;
  maxIterPCG = (iter > maxIterPCG) ? iter : maxIterPCG;

  return info;

}


int JDPCG::solve(int numEigen, Epetra_MultiVector &Q, double *lambda) {

  return JDPCG::reSolve(numEigen, Q, lambda);

}


int JDPCG::reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV) {

  // Computes the smallest eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = M X Lambda
  // 
  // using a Jacobi-Davidson algorithm with PCG (Block version of Notay's algorithm).
  //
  // Reference: "Combination of Jacobi-Davidson and conjugate gradients for the partial
  // symmetric eigenproblem", Y. Notay, Numer. Linear Algebra Appl. (2002), 9:21-44.
  // 
  // Input variables:
  // 
  // numEigen  (integer) = Number of eigenmodes requested
  // 
  // Q (Epetra_MultiVector) = Converged eigenvectors
  //                   The number of columns of Q must be at least numEigen + blockSize.
  //                   The rows of Q are distributed across processors.
  //                   At exit, the first numEigen columns contain the eigenvectors requested.
  // 
  // lambda (array of doubles) = Converged eigenvalues
  //                   At input, it must be of size numEigen + blockSize.
  //                   At exit, the first numEigen locations contain the eigenvalues requested.
  //
  // startingEV (integer) = Number of existing converged eigenvectors
  //                   We assume that the user has check the eigenvectors and
  //                   their M-orthonormality.
  //
  // Return information on status of computation
  // 
  // info >=   0 >> Number of converged eigenpairs at the end of computation
  // 
  // // Failure due to input arguments
  // 
  // info = -  1 >> The stiffness matrix K has not been specified.
  // info = -  2 >> The maps for the matrix K and the matrix M differ.
  // info = -  3 >> The maps for the matrix K and the preconditioner P differ.
  // info = -  4 >> The maps for the vectors and the matrix K differ.
  // info = -  5 >> Q is too small for the number of eigenvalues requested.
  // info = -  6 >> Q is too small for the computation parameters.
  //
  // info = -  7 >> The mass matrix M has not been specified.
  // info = -  8 >> The number of blocks is too small for the number of eigenvalues.
  // 
  // info = - 10 >> Failure during the mass orthonormalization
  // 
  // info = - 30 >> MEMORY
  //

  // Check the input parameters
  
  if (numEigen <= startingEV) {
    return startingEV;
  }

  int info = myVerify.inputArguments(numEigen, K, M, Prec, Q, minimumSpaceDimension(numEigen));
  if (info < 0)
    return info;

  int myPid = MyComm.MyPID();

  if (M == 0) {
    if (myPid == 0) {
      cerr << endl;
      cerr << " !!! The Epetra_Operator object for the mass matrix is not specified !!!" << endl;
      cerr << endl;
    }
    return -7;
  }

  if (numBlock*blockSize < numEigen) {
    if (myPid == 0) {
      cerr << endl;
      cerr << " !!! The space dimension (# of blocks x size of blocks) must be greater than ";
      cerr << " the number of eigenvalues !!!\n";
      cerr << " Number of blocks = " << numBlock << endl;
      cerr << " Size of blocks = " << blockSize << endl;
      cerr << " Number of eigenvalues = " << numEigen << endl;
      cerr << endl;
    }
    return -8;
  }

  // Get the weight for approximating the M-inverse norm
  Epetra_Vector *vectWeight = 0;
  if (normWeight) {
    vectWeight = new Epetra_Vector(View, Q.Map(), normWeight);
  }

  int knownEV = startingEV;
  int localVerbose = verbose*(myPid==0);

  // Define local block vectors
  //
  // PMQ = Preconditioned M-times converged eigenvectors + one working block
  //
  // KX = Working vectors (storing K times one block)
  //
  // R  = Working vectors (storing residuals for one block)
  //
  // MX = Working vectors (storing M times one block)
  //

  int xr = Q.MyLength();
  int dimSearch = blockSize*numBlock;

  Epetra_MultiVector X(View, Q, 0, dimSearch + blockSize);
  if (knownEV > 0) {
    Epetra_MultiVector copyX(View, Q, knownEV, blockSize);
    copyX.Random();
  }
  else {
    X.Random();
  }

  int lwork = 0;
  lwork += (numEigen-startingEV + blockSize)*xr + 2*blockSize*xr;
  lwork += (M) ? blockSize*xr : 0;

  // Workspace for PCG
  lwork += (maxIterLinearSolve > 0) ? 4*xr*blockSize + 6*blockSize + 4*blockSize*blockSize : 0;

  // Workspace for JDCG
  lwork += blockSize + dimSearch + 2*dimSearch*dimSearch;
  lwork += 2*(numEigen-startingEV+blockSize)*(numEigen-startingEV+blockSize);
  lwork += 2*(numEigen-startingEV+blockSize)*blockSize;

  double *work = new (nothrow) double[lwork]; 
  if (work == 0) {
    if (vectWeight)
      delete vectWeight;
    info = -30;
    return info;
  }
  memRequested += sizeof(double)*lwork/(1024.0*1024.0);

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  double *tmpD = work;

  Epetra_MultiVector PMQ(View, Q.Map(), tmpD, xr, numEigen-startingEV+blockSize);
  tmpD = tmpD + (numEigen-startingEV+blockSize)*xr;

  Epetra_MultiVector R(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + blockSize*xr;

  Epetra_MultiVector KX(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + blockSize*xr;

  Epetra_MultiVector MX(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + blockSize*xr;

  // Note: Use MX as workspace in PCG iterations
  double *workPCG = 0;
  if (maxIterLinearSolve > 0) {
    workPCG = tmpD - blockSize*xr;
    tmpD = tmpD + 4*xr*blockSize + 6*blockSize + 4*blockSize*blockSize;
  }

  // theta = Store the local eigenvalues (size: dimSearch)
  //
  // normR = Store the norm of residuals (size: blockSize)
  //
  // KK = Local stiffness matrix         (size: dimSearch x dimSearch)
  //
  // S = Local eigenvectors              (size: dimSearch x dimSearch)
  // 
  // QtMPMQ = Projected "preconditioner" (size: (numEigen-startingEV+blockSize) ^ 2)
  //
  // invQtMPMQ = Inverse of QtMPMQ       (size: (numEigen-startingEV+blockSize) ^ 2)
  //
  // tmpArray = Temporary workspace      (size: 2*(numEigen-startingEV+blockSize) x blockSize)

  double *theta = tmpD;
  tmpD = tmpD + dimSearch;

  double *normR = tmpD;
  tmpD = tmpD + blockSize;

  double *KK = tmpD;
  tmpD = tmpD + dimSearch*dimSearch;
  memset(KK, 0, dimSearch*dimSearch*sizeof(double));

  double *S = tmpD;
  tmpD = tmpD + dimSearch*dimSearch;

  double *QtMPMQ = tmpD;
  tmpD = tmpD + (numEigen-startingEV+blockSize)*(numEigen-startingEV+blockSize);
  int ldQtMPMQ = numEigen - startingEV + blockSize;

  double *invQtMPMQ = tmpD;
  tmpD = tmpD + (numEigen-startingEV+blockSize)*(numEigen-startingEV+blockSize);

  double *tmpArray = tmpD;

  // Define an array to store the residuals history
  if (localVerbose > 2) {
    resHistory = new (nothrow) double[maxIterEigenSolve*blockSize];
    spaceSizeHistory = new (nothrow) int[maxIterEigenSolve];
    iterPCGHistory = new (nothrow) int[maxIterEigenSolve];
    if ((resHistory == 0) || (spaceSizeHistory == 0) || (iterPCGHistory == 0)) {
      if (vectWeight)
        delete vectWeight;
      delete[] work;
      info = -30;
      return info;
    }
    historyCount = 0;
  }

  // Miscellaneous definitions
  bool reStart = false;
  bool criticalExit = false;

  int i, j;
  int nFound = blockSize;
  int bStart = 0;
  int offSet = 0;
  numRestart = 0;

  double tau = 0.0;
  double eta = tau;
  double sqrtTol = sqrt(tolEigenSolve);
  double coefDecay = 0.5;
  double tolPCG = 1.0;

  if (localVerbose > 0) {
    cout << endl;
    cout << " *|* Problem: ";
    if (M)
      cout << "K*Q = M*Q D ";
    else
      cout << "K*Q = Q D ";
    if (Prec)
      cout << " with preconditioner";
    cout << endl;
    cout << " *|* Algorithm = Jacobi-Davidson algorithm with PCG (block version)\n";
    cout << " *|* Size of blocks = " << blockSize << endl;
    cout << " *|* Largest size of search space = " << numBlock*blockSize << endl;
    cout << " *|* Number of requested eigenvalues = " << numEigen << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " *|* Tolerance for convergence = " << tolEigenSolve << endl;
    cout << " *|* Norm used for convergence: ";
    if (vectWeight)
      cout << "weighted L2-norm with user-provided weights" << endl;
    else
      cout << "L^2-norm" << endl;
    if (startingEV > 0)
      cout << " *|* Input converged eigenvectors = " << startingEV << endl;
    cout << "\n -- Start iterations -- \n";
  }

  int maxBlock = (dimSearch/blockSize) - (knownEV/blockSize);

  timeOuterLoop -= MyWatch.WallTime();
  outerIter = 0;
  while (outerIter <= maxIterEigenSolve) {

    highMem = (highMem > currentSize()) ? highMem : currentSize();

    Epetra_MultiVector PMX(View, PMQ, knownEV - startingEV, blockSize);

    int nb;
    for (nb = bStart; nb < maxBlock; ++nb) {

      outerIter += 1;
      if (outerIter > maxIterEigenSolve)
        break;

      int localSize = nb*blockSize;

      Epetra_MultiVector Xcurrent(View, X, localSize + knownEV, blockSize);

      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(Xcurrent, MX);
      timeMassOp += MyWatch.WallTime();
      massOp += blockSize;

      // Orthonormalize X against the known eigenvectors and the previous vectors
      // Note: Use R as a temporary work space
      timeOrtho -= MyWatch.WallTime();
      if (nb == bStart) {
        if (nFound > 0) {
          if (knownEV == 0) {
            info = modalTool.massOrthonormalize(Xcurrent, MX, M, Q, nFound, 2, R.Values());
          }
          else {
            if (localSize == 0) {
              Epetra_MultiVector copyQ(View, X, 0, knownEV);
              info = modalTool.massOrthonormalize(Xcurrent, MX, M, copyQ, nFound, 0, R.Values());
            }
            else {
              Epetra_MultiVector copyQ(View, X, knownEV, localSize);
              info = modalTool.massOrthonormalize(Xcurrent, MX, M, copyQ, nFound, 0, R.Values());
            } // if (localSize == 0)
          } // if (knownEV == 0)
        } // if (nFound > 0)
        nFound = 0;
      }
      else {
        Epetra_MultiVector copyQ(View, X, knownEV, localSize);
        info = modalTool.massOrthonormalize(Xcurrent, MX, M, copyQ, blockSize, 0, R.Values());
      }
      timeOrtho += MyWatch.WallTime();

      // Exit the code when the number of vectors exceeds the space dimension
      if (info < 0) {
        info = -10;
        if (vectWeight)
          delete vectWeight;
        delete[] work;
        if (workPCG)
          delete[] workPCG;
        return info;
      }

      timeStifOp -= MyWatch.WallTime();
      K->Apply(Xcurrent, KX);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;

      if (verbose > 3) {
        if (knownEV + localSize == 0) 
          accuracyCheck(&Xcurrent, &MX, 0);
        else {
          Epetra_MultiVector copyQ(View, X, 0, knownEV + localSize);
          accuracyCheck(&Xcurrent, &MX, &copyQ);
        }
        if (localVerbose > 0)
          cout << endl;
      } // if (verbose > 3)

      // Define the local stiffness matrix
      timeLocalProj -= MyWatch.WallTime();
      for (j = 0; j <= nb; ++j) {
        callBLAS.GEMM('T', 'N', blockSize, blockSize, xr,
                      1.0, X.Values()+(knownEV+j*blockSize)*xr, xr, KX.Values(), xr,
                      0.0, tmpArray + blockSize*blockSize, blockSize);
        MyComm.SumAll(tmpArray + blockSize*blockSize, tmpArray, blockSize*blockSize);
        int iC;
        for (iC = 0; iC < blockSize; ++iC) {
          double *Kpointer = KK + localSize*dimSearch + j*blockSize + iC*dimSearch;
          memcpy(Kpointer, tmpArray + iC*blockSize, blockSize*sizeof(double));
        } // for (iC = 0; iC < blockSize; ++iC)
      } // for (j = 0; j <= nb; ++j)
      timeLocalProj += MyWatch.WallTime();

      // Perform a spectral decomposition
      timeLocalSolve -= MyWatch.WallTime();
      int nevLocal = localSize + blockSize;
      info = modalTool.directSolver(localSize + blockSize, KK, dimSearch, 0, 0,
                                    nevLocal, S, dimSearch, theta, localVerbose, 10);
      timeLocalSolve += MyWatch.WallTime();

      if ((info != 0) && (theta[0] < 0.0)) {
        // Stop as spectral decomposition has a critical failure
        if (info < 0) {
          criticalExit = true;
          break;
        }
        // Restart as spectral decomposition failed
        if (localVerbose > 0) {
          cout << " Iteration " << outerIter;
          cout << "- Failure for spectral decomposition - RESTART with new random search\n";
        }
        reStart = true;
        numRestart += 1;
        timeRestart -= MyWatch.WallTime();
        Epetra_MultiVector Xinit(View, X, knownEV, blockSize);
        Xinit.Random();
        timeRestart += MyWatch.WallTime();
        nFound = blockSize;
        bStart = 0;
        break;
      } // if (info != 0)

      // Start the definition of the new block
      // Note: Use KX as a workspace to store the updated directions
      timeLocalUpdate -= MyWatch.WallTime();
      callBLAS.GEMM('N', 'N', xr, blockSize, localSize+blockSize, 1.0, X.Values()+knownEV*xr, xr,
                    S, dimSearch, 0.0, KX.Values(), xr);
      timeLocalUpdate += MyWatch.WallTime();

      // Apply the mass matrix on the new block
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(KX, MX);
      timeMassOp += MyWatch.WallTime();
      massOp += blockSize;

      // Apply the stiffness matrix on the new block
      timeStifOp -= MyWatch.WallTime();
      K->Apply(KX, R);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;

      // Form the residuals
      timeResidual -= MyWatch.WallTime();
      for (j = 0; j < blockSize; ++j) {
        callBLAS.AXPY(xr, -theta[j], MX.Values()+j*xr, R.Values()+j*xr);
      }
      timeResidual += MyWatch.WallTime();
      residual += blockSize;

      // Compute the norm of residuals
      timeNorm -= MyWatch.WallTime();
      if (vectWeight)
        R.NormWeighted(*vectWeight, normR);
      else
        R.Norm2(normR);
      // Scale the norms of residuals with the eigenvalues
      // Count the number of converged eigenvectors
      nFound = 0;
      for (j = 0; j < blockSize; ++j) {
        normR[j] = (theta[j] == 0.0) ? normR[j] : normR[j]/theta[j];
        if (normR[j] < tolEigenSolve) {
          nFound += 1;
        }
      } // for (j = 0; j < blockSize; ++j)
      timeNorm += MyWatch.WallTime();

      maxSpaceSize = (maxSpaceSize > localSize+blockSize) ? maxSpaceSize : localSize+blockSize;
      sumSpaceSize += localSize + blockSize;

      // Print information on current iteration
      if (localVerbose > 0) {
        cout << " Iteration " << outerIter << " - Number of converged eigenvectors ";
        cout << knownEV + nFound << endl;
      }

      if (localVerbose > 2) {
        memcpy(resHistory + blockSize*historyCount, normR, blockSize*sizeof(double));
        spaceSizeHistory[historyCount] = localSize + blockSize;
      }

      if (localVerbose > 1) {
        cout << endl;
        cout.precision(2);
        cout.setf(ios::scientific, ios::floatfield);
        for (i=0; i<blockSize; ++i) {
          cout << " Iteration " << outerIter << " - Scaled Norm of Residual " << i;
          cout << " = " << normR[i] << endl;
        }
        cout << endl;
        cout.precision(2);
        for (i=0; i<localSize; ++i) {
          cout << " Iteration " << outerIter << " - Ritz eigenvalue " << i;
          cout.setf((fabs(theta[i]) < 0.01) ? ios::scientific : ios::fixed, ios::floatfield);
          cout << " = " << theta[i] << endl;
        }
        cout << endl;
      }

      // Exit the loop when some vectors have converged
      if (nFound > 0) {
        tolPCG = 1.0;
        offSet = 0;
        nb += 1;
        if (localVerbose > 2) {
          iterPCGHistory[historyCount] = 0;
          historyCount += 1;
        }
        break;
      }

      // Exit the loop when all positions are filled
      if ((maxBlock > 1) && (nb == maxBlock - 1)) {
        if (localVerbose > 2) {
          iterPCGHistory[historyCount] = 0;
          historyCount += 1;
        }
        continue;
      }

      // Apply the preconditioner on the new direction
      if (Prec) {
        timePrecOp -= MyWatch.WallTime();
        Prec->ApplyInverse(MX, PMX);
        timePrecOp += MyWatch.WallTime();
        precOp += blockSize;
      } // if (Prec)
      else {
        memcpy(PMX.Values(), MX.Values(), xr*blockSize*sizeof(double));
      } // if (Prec)

      // Update the upper triangular part of the matrix QtMPMQ
      // Note: Use tmpArray as a workspace
      timeBuildQtMPMQ -= MyWatch.WallTime();
      int qLength = knownEV - startingEV + blockSize;
      callBLAS.GEMM('T', 'N', qLength, blockSize, xr, 1.0, PMQ.Values(), xr, MX.Values(), xr,
                    0.0, tmpArray + qLength*blockSize, qLength);
      MyComm.SumAll(tmpArray + qLength*blockSize, tmpArray, qLength*blockSize);
      for (j = 0; j < blockSize; ++j) {
        memcpy(QtMPMQ + (knownEV-startingEV+j)*ldQtMPMQ, tmpArray + j*qLength,
               qLength*sizeof(double));
      }
      timeBuildQtMPMQ += MyWatch.WallTime();

      // Factor the matrix QtMPMQ
      for (j = 0; j < qLength; ++j)
        memcpy(invQtMPMQ + j*ldQtMPMQ, QtMPMQ + j*ldQtMPMQ, (j+1)*sizeof(double));
      callLAPACK.POTRF('U', qLength, invQtMPMQ, ldQtMPMQ, &info);

      // Treat the error messages for Cholesky factorization
      if (info != 0) {
        if (info < 0) {
          if (myPid == 0) {
            cerr << endl;
            cerr << " !!! The argument " << -info << " of DPOTRF had an illegal value !!!\n";
            cerr << endl;
          }
          exit(-1);
        }
        // Restart as factorization failed
        if (localVerbose > 0) {
          cout << " Iteration " << outerIter;
          cout << " - Failure for local factorization (" << info << "/" << qLength << ")";
          cout << " - RESTART with new random search";
          cout << endl;
        }
        reStart = true;
        numRestart += 1;
        timeRestart -= MyWatch.WallTime();
        Epetra_MultiVector Xinit(View, X, knownEV, blockSize);
        Xinit.Random();
        timeRestart += MyWatch.WallTime();
        nFound = blockSize;
        bStart = 0;
        if (localVerbose > 2) {
          iterPCGHistory[historyCount] = 0;
          historyCount += 1;
        }
        break;
      }

      // Correction equation
      // Note that the new working block U is stored in KX,
      // while MX contains M*U and PMX contains P^{-1}*M*U
      Epetra_MultiVector Xnext(View, X, (maxBlock == 1) ? knownEV 
                                                        : knownEV+localSize+blockSize, blockSize);
      if (normR[0] < sqrtTol) 
        eta = theta[0];
      else 
        eta = tau;

      if (verbose > 3) {
        double maxDotQ = 0.0;
        if (knownEV > 0) {
          Epetra_MultiVector copyQ(View, X, 0, knownEV);
          maxDotQ = myVerify.errorOrthogonality(&copyQ, &R, 0);
        }
        double maxDotU = myVerify.errorOrthogonality(&KX, &R, 0);
        if (myPid == 0) {
          cout << " >> Orthogonality for the right hand side of the correction equation = ";
          double tmp = (maxDotU > maxDotQ) ? maxDotU : maxDotQ;
          cout.precision(4);
          cout.setf(ios::scientific, ios::floatfield);
          cout << tmp << endl << endl;
        }
      }

      timeCorrectionSolve -= MyWatch.WallTime();
      if (startingEV > 0) {
        Epetra_MultiVector startQ(View, X, 0, startingEV);
        if (knownEV > startingEV) {
          Epetra_MultiVector copyQ(View, X, startingEV, knownEV - startingEV);
          info = jacobiPCG(R, Xnext, &KX, &copyQ, eta, tolPCG, maxIterLinearSolve, invQtMPMQ,
                           ldQtMPMQ, PMQ.Values(), tmpArray, workPCG, vectWeight, &startQ);
        }
        else {
          info = jacobiPCG(R, Xnext, &KX, 0, eta, tolPCG, maxIterLinearSolve, invQtMPMQ,
                           ldQtMPMQ, PMQ.Values(), tmpArray, workPCG, vectWeight, &startQ);
        }
      }
      else {
        if (knownEV > 0) {
          Epetra_MultiVector copyQ(View, X, 0, knownEV);
          info = jacobiPCG(R, Xnext, &KX, &copyQ, eta, tolPCG, maxIterLinearSolve, invQtMPMQ,
                           ldQtMPMQ, PMQ.Values(), tmpArray, workPCG, vectWeight);
        }
        else {
          info = jacobiPCG(R, Xnext, &KX, 0, eta, tolPCG, maxIterLinearSolve, invQtMPMQ,
                           ldQtMPMQ, PMQ.Values(), tmpArray, workPCG, vectWeight);
        }
      }
      timeCorrectionSolve += MyWatch.WallTime();

      if (verbose > 3) {
        double maxDotQ = 0.0;
        if (knownEV > 0) {
          Epetra_MultiVector copyQ(View, X, 0, knownEV);
          maxDotQ = myVerify.errorOrthogonality(&copyQ, &Xnext, M);
        }
        double maxDotU = myVerify.errorOrthogonality(&KX, &Xnext, M);
        if (myPid == 0) {
          double tmp = (maxDotU > maxDotQ) ? maxDotU : maxDotQ;
          cout << " >> Orthogonality for the solution of the correction equation = ";
          cout.precision(4);
          cout.setf(ios::scientific, ios::floatfield);
          cout << tmp << endl << endl;
        }
      }

      numCorrectionSolve += blockSize;
      if (info < 0) {
        if ((info == -1) || (info == -maxIterLinearSolve-1)) {
          if (localVerbose > 0) {
            cout << " Iteration " << outerIter;
            cout << " - Failure inside PCG";
            cout << " - RESTART with new random search";
            cout << endl;
          }
          if (localVerbose > 2) {
            iterPCGHistory[historyCount] = -1;
            historyCount += 1;
          }
          reStart = true;
          numRestart += 1;
          timeRestart -= MyWatch.WallTime();
          Epetra_MultiVector Xinit(View, X, knownEV, blockSize);
          Xinit.Random();
          timeRestart += MyWatch.WallTime();
          nFound = blockSize;
          bStart = 0;
          info = 0;
          break;
        } // if ((info == -1) || (info == -iterMax-1))
        else {
          if (localVerbose > 2) {
            iterPCGHistory[historyCount] = (info < 0) ? -info : info;
            historyCount += 1;
          }
        } // if ((info == -1) || (info == -iterMax-1))
      } // if (info < 0)
      else {
        if (localVerbose > 2) {
          iterPCGHistory[historyCount] = info;
          historyCount += 1;
        }
      } // if (info < 0)
      info = 0;

      tolPCG *= coefDecay;

      if (maxBlock == 1) {
        Xcurrent.Update(1.0, KX, -1.0);
      }

    } // for (nb = bstart; nb < maxBlock; ++nb)

    if (outerIter > maxIterEigenSolve)
      break;

    if (reStart == true) {
      reStart = false;
      continue;
    }

    if (criticalExit == true)
      break;

    // Store the final converged eigenvectors
    if (knownEV + nFound >= numEigen) {
      for (j = 0; j < blockSize; ++j) {
        if (normR[j] < tolEigenSolve) {
          memcpy(X.Values() + knownEV*xr, KX.Values() + j*xr, xr*sizeof(double));
          lambda[knownEV] = theta[j];
          knownEV += 1;
        }
      }
      if (localVerbose == 1) {
        cout << endl;
        cout.precision(2);
        cout.setf(ios::scientific, ios::floatfield);
        for (i=0; i<blockSize; ++i) {
          cout << " Iteration " << outerIter << " - Scaled Norm of Residual " << i;
          cout << " = " << normR[i] << endl;
        }
        cout << endl;
      }  
      break;
    } // if (knownEV + nFound >= numEigen)

    // Treat the particular case of one block
    if ((maxBlock == 1) && (nFound == 0)) {
      nFound = blockSize;
      continue;
    }

    // Define the restarting space when there is only one block
    if (maxBlock == 1) {
      timeRestart -= MyWatch.WallTime();
      double *Xpointer = X.Values() + (knownEV+nFound)*xr;
      nFound = 0;
      for (j = 0; j < blockSize; ++j) {
        if (normR[j] < tolEigenSolve) {
          memcpy(X.Values() + knownEV*xr, KX.Values() + j*xr, xr*sizeof(double));
          lambda[knownEV] = theta[j];
          knownEV += 1;
          nFound += 1;
        }
        else {
          memcpy(Xpointer + (j-nFound)*xr, KX.Values() + j*xr, xr*sizeof(double));
        }
      }
      knownEV -= nFound;
      Epetra_MultiVector Xnext(View, X, knownEV + blockSize, nFound);
      Xnext.Random();
      timeRestart += MyWatch.WallTime();
    }

    // Define the restarting block when there is more than one block
    int oldCol, newCol;
    if (maxBlock > 1) {
      timeRestart -= MyWatch.WallTime();
      int firstIndex = blockSize;
      for (j = 0; j < blockSize; ++j) {
        if (normR[j] >= tolEigenSolve) {
          firstIndex = j;
          break;
        }
      } // for (j = 0; j < blockSize; ++j)
      while (firstIndex < nFound) {
        for (j = firstIndex; j < blockSize; ++j) {
          if (normR[j] < tolEigenSolve) {
            // Swap the j-th and firstIndex-th position
            callFortran.SWAP(nb*blockSize, S + j*dimSearch, 1, S + firstIndex*dimSearch, 1);
            callFortran.SWAP(1, theta + j, 1, theta + firstIndex, 1);
            callFortran.SWAP(1, normR + j, 1, normR + firstIndex, 1);
            break;
          }
        } // for (j = firstIndex; j < blockSize; ++j)
        for (j = 0; j < blockSize; ++j) {
          if (normR[j] >= tolEigenSolve) {
            firstIndex = j;
            break;
          }
        } // for (j = 0; j < blockSize; ++j)
      } // while (firstIndex < nFound)

      // Copy the converged eigenvalues
      memcpy(lambda + knownEV, theta, nFound*sizeof(double));

      // Define the restarting size
      bStart = ((nb - offSet) > 2) ? (nb - offSet)/2 : 0;

      // Define the restarting space and local stiffness
      memset(KK, 0, nb*blockSize*dimSearch*sizeof(double));
      for (j = 0; j < bStart*blockSize; ++j) {
        KK[j + j*dimSearch] = theta[j + nFound];
      }
      // Form the restarting space
      oldCol = nb*blockSize;
      newCol = nFound + (bStart+1)*blockSize;
      newCol = (newCol > oldCol) ? oldCol : newCol;
      callFortran.GEQRF(oldCol, newCol, S, dimSearch, theta, R.Values(), xr*blockSize, &info);
      callFortran.ORMQR('R', 'N', xr, oldCol, newCol, S, dimSearch, theta,
                        X.Values()+knownEV*xr, xr, R.Values(), xr*blockSize, &info);
      // Put random vectors if the Rayleigh Ritz vectors are not enough
      newCol = nFound + (bStart+1)*blockSize;
      if (newCol > oldCol) {
        Epetra_MultiVector Xnext(View, X, knownEV+blockSize-nFound, nFound);
        Xnext.Random();
      }
      timeRestart += MyWatch.WallTime();
    } // if (maxBlock > 1)

    if (nFound > 0) {
      // Update MQ
      Epetra_MultiVector Qconv(View, X, knownEV, nFound);
      Epetra_MultiVector MQconv(View, MX, 0, nFound);
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(Qconv, MQconv);
      timeMassOp += MyWatch.WallTime();
      massOp += nFound;
      // Update PMQ
      Epetra_MultiVector PMQconv(View, PMQ, knownEV-startingEV, nFound);
      if (Prec) {
        timePrecOp -= MyWatch.WallTime();
        Prec->ApplyInverse(MQconv, PMQconv);
        timePrecOp += MyWatch.WallTime();
        precOp += nFound;
      }
      else {
        memcpy(PMQconv.Values(), MQconv.Values(), xr*nFound*sizeof(double));
      }
      // Update QtMPMQ
      // Note: Use tmpArray as workspace
      int newEV = knownEV - startingEV;
      timeBuildQtMPMQ -= MyWatch.WallTime();
      callBLAS.GEMM('T', 'N', newEV + nFound, nFound, xr, 1.0, PMQ.Values(), xr,
                    MX.Values(), xr, 0.0, QtMPMQ + newEV*ldQtMPMQ, newEV + nFound);
      MyComm.SumAll(QtMPMQ + newEV*ldQtMPMQ, tmpArray, (newEV + nFound)*nFound);
      for (j = 0; j < nFound; ++j) {
        memcpy(QtMPMQ + (newEV+j)*ldQtMPMQ, tmpArray + j*(newEV+nFound), 
               (newEV+nFound)*sizeof(double));
      }
      timeBuildQtMPMQ += MyWatch.WallTime();
    } // if (nFound > 0)
    else {
      offSet += 1;
    } // if (nFound > 0)

    knownEV += nFound;
    maxBlock = (dimSearch/blockSize) - (knownEV/blockSize);

    // The value of nFound commands how many vectors will be orthogonalized.
    if ((maxBlock > 1) && (newCol <= oldCol))
      nFound = 0;

  } // while (outerIter <= maxIterEigenSolve)
  timeOuterLoop += MyWatch.WallTime();
  highMem = (highMem > currentSize()) ? highMem : currentSize();

  // Clean memory
  delete[] work;
  if (vectWeight)
    delete vectWeight;

  // Sort the eigenpairs
  timePostProce -= MyWatch.WallTime();
  if ((info == 0) && (knownEV > 0)) {
    mySort.sortScalars_Vectors(knownEV, lambda, Q.Values(), Q.MyLength());
  }
  timePostProce += MyWatch.WallTime();

  return (info == 0) ? knownEV : info;

}


void JDPCG::accuracyCheck(const Epetra_MultiVector *X, const Epetra_MultiVector *MX,
                             const Epetra_MultiVector *Q) const {

  cout.precision(2);
  cout.setf(ios::scientific, ios::floatfield);
  double tmp;

  int myPid = MyComm.MyPID();

  if (X) {
    if (M) {
      if (MX) {
        tmp = myVerify.errorEquality(X, MX, M);
        if (myPid == 0)
          cout << " >> Difference between MX and M*X = " << tmp << endl;
      }
      tmp = myVerify.errorOrthonormality(X, M);
      if (myPid == 0)
        cout << " >> Error in X^T M X - I = " << tmp << endl;
    }
    else {
      tmp = myVerify.errorOrthonormality(X, 0);
      if (myPid == 0)
        cout << " >> Error in X^T X - I = " << tmp << endl;
    }
  }

  if (Q == 0)
    return;

  if (M) {
    tmp = myVerify.errorOrthonormality(Q, M);
    if (myPid == 0)
      cout << " >> Error in Q^T M Q - I = " << tmp << endl;
    if (X) {
      tmp = myVerify.errorOrthogonality(Q, X, M);
      if (myPid == 0)
        cout << " >> Orthogonality Q^T M X up to " << tmp << endl;
    }
  }
  else {
    tmp = myVerify.errorOrthonormality(Q, 0);
    if (myPid == 0)
      cout << " >> Error in Q^T Q - I = " << tmp << endl;
    if (X) {
      tmp = myVerify.errorOrthogonality(Q, X, 0);
      if (myPid == 0)
        cout << " >> Orthogonality Q^T X up to " << tmp << endl;
    }
  }

}


int JDPCG::minimumSpaceDimension(int nev) const {

  int myPid = MyComm.MyPID();
  
  if ((myPid == 0) && (numBlock*blockSize < nev)) {
    cerr << endl;
    cerr << " !!! The space dimension (# of blocks x size of blocks) must be greater than ";
    cerr << " the number of eigenvalues !!!\n";
    cerr << " Number of blocks = " << numBlock << endl;
    cerr << " Size of blocks = " << blockSize << endl;
    cerr << " Number of eigenvalues = " << nev << endl;
    cerr << endl;
    return -1;
  }

  return nev + blockSize;

}


void JDPCG::initializeCounters() {

  historyCount = 0;
  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

  maxSpaceSize = 0;
  sumSpaceSize = 0;
  if (spaceSizeHistory) {
    delete[] spaceSizeHistory;
    spaceSizeHistory = 0;
  }

  maxIterPCG = 0;
  sumIterPCG = 0;
  if (iterPCGHistory) {
    delete[] iterPCGHistory;
    iterPCGHistory = 0;
  }

  memRequested = 0.0;
  highMem = 0.0;

  massOp = 0;
  numCorrectionPrec = 0;
  numCorrectionSolve = 0;
  numPCGmassOp = 0;
  numPCGstifOp = 0;
  numRestart = 0;
  outerIter = 0;
  precOp = 0;
  residual = 0;
  stifOp = 0;

  timeBuildQtMPMQ = 0.0;
  timeCorrectionPrec = 0.0;
  timeCorrectionSolve = 0.0;
  timeLocalProj = 0.0;
  timeLocalSolve = 0.0;
  timeLocalUpdate = 0.0;
  timeMassOp = 0.0;
  timeNorm = 0.0;
  timeOrtho = 0.0;
  timeOuterLoop = 0.0;
  timePCGEigCheck = 0.0;
  timePCGLoop = 0.0;
  timePCGOpMult = 0.0;
  timePCGPrec = 0.0;
  timePostProce = 0.0;
  timePrecOp = 0.0;
  timeResidual = 0.0;
  timeRestart = 0.0;
  timeStifOp = 0.0;

}


void JDPCG::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: Jacobi-Davidson algorithm with PCG (block version)\n";
    cout << " Block Size: " << blockSize << endl;
    cout << " Number of Blocks kept: " << numBlock << endl;
  }

}


void JDPCG::historyInfo() const {

  if (resHistory) {
    int j;
    cout << " Block Size: " << blockSize << endl;
    cout << endl;
    cout << " Residuals    Search Space Dim.   Inner Iter. \n";
    cout << endl;
    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);
    for (j = 0; j < historyCount; ++j) {
      int ii;
      for (ii = 0; ii < blockSize; ++ii) {
        cout << resHistory[blockSize*j + ii] << "      ";
        cout.width(4);
        cout << spaceSizeHistory[j] << "          ";
        cout.width(4);
        cout << iterPCGHistory[j] << endl;
      }
    }
    cout << endl;
  }

}


void JDPCG::memoryInfo() const {

  int myPid = MyComm.MyPID();

  double maxHighMem = 0.0;
  double tmp = highMem;
  MyComm.MaxAll(&tmp, &maxHighMem, 1);

  double maxMemRequested = 0.0;
  tmp = memRequested;
  MyComm.MaxAll(&tmp, &maxMemRequested, 1);

  if (myPid == 0) {
    cout.precision(2);
    cout.setf(ios::fixed, ios::floatfield);
    cout << " Memory requested per processor by the eigensolver   = (EST) ";
    cout.width(6);
    cout << maxMemRequested << " MB " << endl;
    cout << endl;
    cout << " High water mark in eigensolver                      = (EST) ";
    cout.width(6);
    cout << maxHighMem << " MB " << endl;
    cout << endl;
  }

}


void JDPCG::operationInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Operations ---\n\n";
    cout << " Total number of mass matrix multiplications      = ";
    cout.width(9);
    int tmp = massOp + modalTool.getNumProj_MassMult() + modalTool.getNumNorm_MassMult();
    tmp += numPCGmassOp;
    cout << tmp << endl;
    cout << "       Mass multiplications in solver loop        = ";
    cout.width(9);
    cout << massOp << endl;
    cout << "       Mass multiplications for orthogonalization = ";
    cout.width(9);
    cout << modalTool.getNumProj_MassMult() << endl;
    cout << "       Mass multiplications for normalization     = ";
    cout.width(9);
    cout << modalTool.getNumNorm_MassMult() << endl;
    cout << "       Mass multiplications in PCG loop           = ";
    cout.width(9);
    cout << numPCGmassOp << endl;
    cout << endl;
    cout << " Total number of stiffness matrix operations      = ";
    cout.width(9);
    cout << stifOp + numPCGstifOp << endl;
    cout << "       Stiffness multiplications in solver loop   = ";
    cout.width(9);
    cout << stifOp << endl;
    cout << "       Stiffness multiplications in PCG loop      = ";
    cout.width(9);
    cout << numPCGstifOp << endl;
    cout << endl;
    cout << " Total number of preconditioner applications      = ";
    cout.width(9);
    cout << precOp << endl;
    cout << endl;
    cout << " Total number of PCG solve                        = ";
    cout.width(9);
    cout << numCorrectionSolve << endl;
    cout << "       Number of projected precond. appl.  : ";
    cout.width(9);
    cout << numCorrectionPrec << endl;
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(2);
    cout.width(9);
    cout << "       Average number of iter. per solve   : ";
    cout.width(9);
    cout << ((double) sumIterPCG)*blockSize/((double) numCorrectionSolve) << endl;
    cout << "       Maximum number of iter. per solve   : ";
    cout.width(9);
    cout << maxIterPCG << endl;
    cout << endl;
    cout << " Total number of computed eigen-residuals         = ";
    cout.width(9);
    cout << residual << endl;
    cout << endl;
    cout << " Total number of outer iterations                 = ";
    cout.width(9);
    cout << outerIter << endl;
    cout << "       Number of restarts                         = ";
    cout.width(9);
    cout << numRestart << endl;
    cout << endl;
    cout << " Maximum size of search space                     = ";
    cout.width(9);
    cout << maxSpaceSize << endl;
    cout << " Average size of search space                     = ";
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(2);
    cout.width(9);
    cout << ((double) sumSpaceSize)/((double) outerIter) << endl;
    cout << endl;
  }

}


void JDPCG::timeInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Timings ---\n\n";
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(2);
    cout << " Total time for outer loop                       = ";
    cout.width(9);
    cout << timeOuterLoop << " s" << endl;
    cout << "       Time for mass matrix operations           = ";
    cout.width(9);
    cout << timeMassOp << " s     ";
    cout.width(5);
    cout << 100*timeMassOp/timeOuterLoop << " %\n";
    cout << "       Time for stiffness matrix operations      = ";
    cout.width(9);
    cout << timeStifOp << " s     ";
    cout.width(5);
    cout << 100*timeStifOp/timeOuterLoop << " %\n";
    cout << "       Time for preconditioner applications      = ";
    cout.width(9);
    cout << timePrecOp << " s     ";
    cout.width(5);
    cout << 100*timePrecOp/timeOuterLoop << " %\n";
    cout << "       Time for orthogonalizations               = ";
    cout.width(9);
    cout << timeOrtho << " s     ";
    cout.width(5);
    cout << 100*timeOrtho/timeOuterLoop << " %\n";
    cout << "            Projection step          : ";
    cout.width(9);
    cout << modalTool.getTimeProj() << " s\n";
    cout << "            Normalization step       : ";
    cout.width(9);
    cout << modalTool.getTimeNorm() << " s\n";
    cout << "       Time for local projection                 = ";
    cout.width(9);
    cout << timeLocalProj << " s     ";
    cout.width(5);
    cout << 100*timeLocalProj/timeOuterLoop << " %\n";
    cout << "       Time for local eigensolve                 = ";
    cout.width(9);
    cout << timeLocalSolve << " s     ";
    cout.width(5);
    cout << 100*timeLocalSolve/timeOuterLoop << " %\n";
    cout << "       Time for local update                     = ";
    cout.width(9);
    cout << timeLocalUpdate << " s     ";
    cout.width(5);
    cout << 100*timeLocalUpdate/timeOuterLoop << " %\n";
    cout << "       Time for building the matrix QtMP^{-1}MQ  = ";
    cout.width(9);
    cout << timeBuildQtMPMQ << " s     ";
    cout.width(5);
    cout << 100*timeBuildQtMPMQ/timeOuterLoop << " %\n";
    cout << "       Time for solving the correction equation  = ";
    cout.width(9);
    cout << timeCorrectionSolve << " s     ";
    cout.width(5);
    cout << 100*timeCorrectionSolve/timeOuterLoop << " %\n";
    cout << "            Mult. Preconditioner     : ";
    cout.width(9);
    cout << timeCorrectionPrec << " s" << endl;
    cout << "            Correction of Prec.      : ";
    cout.width(9);
    cout << (timePCGPrec - timeCorrectionPrec) << " s" << endl;
    cout << "            Shifted Matrix Mult.     : ";
    cout.width(9);
    cout << timePCGOpMult << " s" << endl;
    cout << "            Eigen-residuals checks   : ";
    cout.width(9);
    cout << timePCGEigCheck << " s" << endl;
    cout << "       Time for restarting space definition      = ";
    cout.width(9);
    cout << timeRestart << " s     ";
    cout.width(5);
    cout << 100*timeRestart/timeOuterLoop << " %\n";
    cout << "       Time for residual computations            = ";
    cout.width(9);
    cout << timeResidual << " s     ";
    cout.width(5);
    cout << 100*timeResidual/timeOuterLoop << " %\n";
    cout << "\n";
    cout << " Total time for post-processing                  = ";
    cout.width(9);
    cout << timePostProce << " s\n";
    cout << endl;
  } // if (myPid == 0)

}


