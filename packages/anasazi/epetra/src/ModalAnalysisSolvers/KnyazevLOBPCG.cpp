// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
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
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#include "KnyazevLOBPCG.h"


KnyazevLOBPCG::KnyazevLOBPCG(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                                   const Epetra_Operator *PP,
                                   double _tol, int _maxIter, int _verb) :
           myVerify(_Comm),
           callBLAS(),
           callFortran(),
           modalTool(_Comm),
           mySort(),
           MyComm(_Comm),
           K(KK),
           M(0),
           Prec(PP),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           blockSize(0),
           normWeight(0),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           numRestart(0),
           outerIter(0),
           precOp(0),
           residual(0),
           stifOp(0),
           timeLocalProj(0.0),
           timeLocalSolve(0.0),
           timeLocalUpdate(0.0),
           timeMassOp(0.0),
           timeNorm(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timePrecOp(0.0),
           timeResidual(0.0),
           timeStifOp(0.0)
           {

}


KnyazevLOBPCG::KnyazevLOBPCG(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                                   const Epetra_Operator *MM, const Epetra_Operator *PP,
                                   double _tol, int _maxIter, int _verb,
                                   double *_weight) :
           myVerify(_Comm),
           callBLAS(),
           callFortran(),
           modalTool(_Comm),
           mySort(),
           MyComm(_Comm),
           K(KK),
           M(MM),
           Prec(PP),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           blockSize(0),
           normWeight(_weight),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           numRestart(0),
           outerIter(0),
           precOp(0),
           residual(0),
           stifOp(0),
           timeLocalProj(0.0),
           timeLocalSolve(0.0),
           timeLocalUpdate(0.0),
           timeMassOp(0.0),
           timeNorm(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timePrecOp(0.0),
           timeResidual(0.0),
           timeStifOp(0.0)
           {

}


KnyazevLOBPCG::~KnyazevLOBPCG() {

  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

}


int KnyazevLOBPCG::solve(int numEigen, Epetra_MultiVector &Q, double *lambda) {

  return KnyazevLOBPCG::reSolve(numEigen, Q, lambda);

}


int KnyazevLOBPCG::reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV) {

  // Computes the smallest eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = M X Lambda
  // 
  // using a Locally Optimal Block Preconditioned Conjugate Gradient method. 
  //
  // Note that if M is not specified, then  K X = X Lambda is solved.
  // Note that the blocksize is equal to the number of requested eigenvalues.
  // 
  // Ref: A. Knyazev, "Toward the optimal preconditioned eigensolver:
  // Locally optimal block preconditioner conjugate gradient method",
  // SIAM J. Sci. Comput., vol 23, n 2, pp. 517-541
  // Ref: A. Knyazev and M. Argentati, "Implementation of a preconditioned
  // eigensolver using Hypre", Numerical Linear Algebra with Applications (submitted)
  // 
  // Input variables:
  // 
  // numEigen  (integer) = Number of eigenmodes requested
  // 
  // Q (Epetra_MultiVector) = Converged eigenvectors
  //                   The number of columns of Q must be equal to numEigen.
  //                   The rows of Q are distributed across processors.
  //                   At exit, the first numEigen columns contain the eigenvectors requested.
  // 
  // lambda (array of doubles) = Converged eigenvalues
  //                   At input, it must be of size numEigen.
  //                   At exit, the first numEigen locations contain the eigenvalues requested.
  //
  // startingEV (integer) = Number of existing converged eigenmodes
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
  // info = - 10 >> Failure during the mass orthonormalization
  // 
  // info = - 20 >> Error in LAPACK during the local eigensolve
  //
  // info = - 30 >> MEMORY
  //

  // Check the input parameters
  
  if (numEigen <= startingEV) {
    return startingEV;
  }

  int info = myVerify.inputArguments(numEigen, K, M, Prec, Q, numEigen);
  if (info < 0)
    return info;

  int myPid = MyComm.MyPID();

  // Get the weight for approximating the M-inverse norm
  Epetra_Vector *vectWeight = 0;
  if (normWeight) {
    vectWeight = new Epetra_Vector(View, Q.Map(), normWeight);
  }

  int knownEV = startingEV;
  int localVerbose = verbose*(myPid==0);

  // Define local block vectors
  //
  // MX = Working vectors (storing M*X if M is specified, else pointing to X)
  // KX = Working vectors (storing K*X)
  //
  // R = Residuals
  //
  // H = Preconditioned search space
  // MH = Working vectors (storing M*H if M is specified, else pointing to H)
  // KH = Working vectors (storing K*H)
  //
  // P = Search directions
  // MP = Working vectors (storing M*P if M is specified, else pointing to P)
  // KP = Working vectors (storing K*P)

  int xr = Q.MyLength();
  Epetra_MultiVector X(View, Q, 0, numEigen);

  blockSize = numEigen;

  int tmp;
  tmp = (M == 0) ? 6*numEigen*xr : 9*numEigen*xr;

  double *work1 = new (nothrow) double[tmp]; 
  if (work1 == 0) {
    if (vectWeight)
      delete vectWeight;
    info = -30;
    return info;
  }
  memRequested += sizeof(double)*tmp/(1024.0*1024.0);

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  double *tmpD = work1;

  Epetra_MultiVector KX(View, Q.Map(), tmpD, xr, numEigen);
  tmpD = tmpD + xr*numEigen;

  Epetra_MultiVector MX(View, Q.Map(), (M) ? tmpD : X.Values(), xr, numEigen);
  tmpD = (M) ? tmpD + xr*numEigen : tmpD;

  Epetra_MultiVector R(View, Q.Map(), tmpD, xr, numEigen);
  tmpD = tmpD + xr*numEigen;

  Epetra_MultiVector H(View, Q.Map(), tmpD, xr, numEigen);
  tmpD = tmpD + xr*numEigen;

  Epetra_MultiVector KH(View, Q.Map(), tmpD, xr, numEigen);
  tmpD = tmpD + xr*numEigen;

  Epetra_MultiVector MH(View, Q.Map(), (M) ? tmpD : H.Values(), xr, numEigen);
  tmpD = (M) ? tmpD + xr*numEigen : tmpD;

  Epetra_MultiVector P(View, Q.Map(), tmpD, xr, numEigen);
  tmpD = tmpD + xr*numEigen;

  Epetra_MultiVector KP(View, Q.Map(), tmpD, xr, numEigen);
  tmpD = tmpD + xr*numEigen;

  Epetra_MultiVector MP(View, Q.Map(), (M) ? tmpD : P.Values(), xr, numEigen);

  if (startingEV > 0) {
    // Fill the first vectors of KX and MX
    Epetra_MultiVector copyX(View, X, 0, startingEV);
    Epetra_MultiVector copyKX(View, KX, 0, startingEV);
    timeStifOp -= MyWatch.WallTime();
    K->Apply(copyX, copyKX); 
    timeStifOp += MyWatch.WallTime();
    stifOp += startingEV;
    timeMassOp -= MyWatch.WallTime();
    if (M) {
      Epetra_MultiVector copyMX(View, MX, 0, startingEV);
      M->Apply(copyX, copyMX); 
    }
    timeMassOp += MyWatch.WallTime();
    massOp += startingEV;
  }

  // Define arrays
  //
  // theta = Store the local eigenvalues (size: numEigen)
  // normR = Store the norm of residuals (size: numEigen)
  //
  // MM = Local mass matrix              (size: 3*numEigen x 3*numEigen)
  // KK = Local stiffness matrix         (size: 3*numEigen x 3*numEigen)
  //
  // S = Local eigenvectors              (size: 3*numEigen x 3*numEigen)

  int lwork2;
  lwork2 = 2*numEigen + 27*numEigen*numEigen;
  double *work2 = new (nothrow) double[lwork2];
  if (work2 == 0) {
    if (vectWeight)
      delete vectWeight;
    delete[] work1;
    info = -30;
    return info;
  }

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  tmpD = work2;

  double *theta = tmpD;
  tmpD = tmpD + numEigen;

  double *normR = tmpD;
  tmpD = tmpD + numEigen;

  double *MM = tmpD;
  tmpD = tmpD + 9*numEigen*numEigen;

  double *KK = tmpD;
  tmpD = tmpD + 9*numEigen*numEigen;

  double *S = tmpD;

  memRequested += sizeof(double)*lwork2/(1024.0*1024.0);

  // Define an array to store the residuals history
  if (localVerbose > 2) {
    resHistory = new (nothrow) double[maxIterEigenSolve*numEigen];
    if (resHistory == 0) {
      if (vectWeight)
        delete vectWeight;
      delete[] work1;
      delete[] work2;
      info = -30;
      return info;
    }
    historyCount = 0;
  }

  // Miscellaneous definitions

  bool reStart = false;
  numRestart = 0;

  int localSize;
  int firstIndex = knownEV;
  int i, j;

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
    cout << " *|* Algorithm = LOBPCG (Knyazev's version)" << endl;
    cout << " *|* Size of blocks = " << blockSize << endl;
    cout << " *|* Number of requested eigenvalues = " << numEigen << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " *|* Tolerance for convergence = " << tolEigenSolve << endl;
    cout << " *|* Norm used for convergence: ";
    if (normWeight)
      cout << "weighted L2-norm with user-provided weights" << endl;
    else
      cout << "L^2-norm" << endl;
    if (startingEV > 0)
      cout << " *|* Input converged eigenvectors = " << startingEV << endl;
    cout << "\n -- Start iterations -- \n";
  }

  timeOuterLoop -= MyWatch.WallTime();
  for (outerIter = 1; outerIter <= maxIterEigenSolve; ++outerIter) {

    highMem = (highMem > currentSize()) ? highMem : currentSize();

    int workingSize = numEigen - knownEV;

    Epetra_MultiVector  XI(View, X, firstIndex, workingSize);
    Epetra_MultiVector MXI(View, MX, firstIndex, workingSize);
    Epetra_MultiVector KXI(View, KX, firstIndex, workingSize);

    Epetra_MultiVector  HI(View, H, firstIndex, workingSize);
    Epetra_MultiVector MHI(View, MH, firstIndex, workingSize);
    Epetra_MultiVector KHI(View, KH, firstIndex, workingSize);

    Epetra_MultiVector  PI(View, P, firstIndex, workingSize);
    Epetra_MultiVector MPI(View, MP, firstIndex, workingSize);
    Epetra_MultiVector KPI(View, KP, firstIndex, workingSize);

    Epetra_MultiVector  RI(View, R, firstIndex, workingSize);

    if ((outerIter == 1) || (reStart == true)) {

      reStart = false;
      localSize = numEigen;

      // Apply the mass matrix to X
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(XI, MXI);
      timeMassOp += MyWatch.WallTime();
      massOp += workingSize;

      // Apply the stiffness matrix to X
      timeStifOp -= MyWatch.WallTime();
      K->Apply(XI, KXI);
      timeStifOp += MyWatch.WallTime();
      stifOp += workingSize;

    } // if ((outerIter == 1) || (reStart == true))
    else {

      // Apply the preconditioner on the residuals
      if (Prec) {
        timePrecOp -= MyWatch.WallTime();
        Prec->ApplyInverse(RI, HI);
        timePrecOp += MyWatch.WallTime();
        precOp += workingSize;
      }
      else {
        memcpy(HI.Values(), RI.Values(), xr*workingSize*sizeof(double));
      }

      // Apply the mass matrix on H
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(HI, MHI);
      timeMassOp += MyWatch.WallTime();
      massOp += workingSize;

      // Apply the stiffness matrix to H
      timeStifOp -= MyWatch.WallTime();
      K->Apply(HI, KHI);
      timeStifOp += MyWatch.WallTime();
      stifOp += workingSize;

      if (localSize == numEigen)
        localSize += workingSize;

    } // if ((outerIter == 1) || (reStart==true))

    // Form "local" mass and stiffness matrices
    // Note: Use S as a temporary workspace
    timeLocalProj -= MyWatch.WallTime();
    modalTool.localProjection(numEigen, numEigen, xr, X.Values(), xr, KX.Values(), xr, 
                              KK, localSize, S);
    modalTool.localProjection(numEigen, numEigen, xr, X.Values(), xr, MX.Values(), xr, 
                              MM, localSize, S);
    if (localSize > numEigen) {
      double *pointer = KK + numEigen*localSize;
      modalTool.localProjection(numEigen, workingSize, xr, X.Values(), xr, KHI.Values(), xr,
                                pointer, localSize, S);
      modalTool.localProjection(workingSize, workingSize, xr, HI.Values(), xr, KHI.Values(), xr,
                                pointer + numEigen, localSize, S);
      pointer = MM + numEigen*localSize;
      modalTool.localProjection(numEigen, workingSize, xr, X.Values(), xr, MHI.Values(), xr, 
                                pointer, localSize, S);
      modalTool.localProjection(workingSize, workingSize, xr, HI.Values(), xr, MHI.Values(), xr,
                                pointer + numEigen, localSize, S);
      if (localSize > numEigen + workingSize) {
        pointer = KK + (numEigen + workingSize)*localSize;
        modalTool.localProjection(numEigen, workingSize, xr, X.Values(), xr, KPI.Values(), xr,
                                  pointer, localSize, S);
        modalTool.localProjection(workingSize, workingSize, xr, HI.Values(), xr, KPI.Values(), xr,
                                  pointer + numEigen, localSize, S);
        modalTool.localProjection(workingSize, workingSize, xr, PI.Values(), xr, KPI.Values(), xr,
                                  pointer + numEigen + workingSize, localSize, S);
        pointer = MM + (numEigen + workingSize)*localSize;
        modalTool.localProjection(numEigen, workingSize, xr, X.Values(), xr, MPI.Values(), xr,
                                  pointer, localSize, S);
        modalTool.localProjection(workingSize, workingSize, xr, HI.Values(), xr, MPI.Values(), xr,
                                  pointer + numEigen, localSize, S);
        modalTool.localProjection(workingSize, workingSize, xr, PI.Values(), xr, MPI.Values(), xr,
                                  pointer + numEigen + workingSize, localSize, S);
      } // if (localSize > numEigen + workingSize)
    } // if (localSize > numEigen)
    timeLocalProj += MyWatch.WallTime();

    // Perform a spectral decomposition
    timeLocalSolve -= MyWatch.WallTime();
    int nevLocal = localSize;
    info = modalTool.directSolver(localSize, KK, localSize, MM, localSize, nevLocal,
                                  S, localSize, theta, localVerbose, 
                                  (blockSize == 1) ? 1 : 0);
    timeLocalSolve += MyWatch.WallTime();

    if (info < 0) {
      // Stop when spectral decomposition has a critical failure
      break;
    } // if (info < 0)

    // Check for restarting
    if ((theta[0] < 0.0) || (nevLocal < numEigen)) {
      if (localVerbose > 0) {
        cout << " Iteration " << outerIter;
        cout << " - Failure for spectral decomposition - RESTART with new random search\n";
      }
      if (workingSize == 1) {
        XI.Random();
      }
      else {
        Epetra_MultiVector Xinit(View, XI, 1, workingSize-1);
        Xinit.Random();
      } // if (workingSize == 1)
      reStart = true;
      numRestart += 1;
      info = 0;
      continue;
    } // if ((theta[0] < 0.0) || (nevLocal < numEigen))

    if ((localSize == numEigen+workingSize) && (nevLocal == numEigen)) {
      for (j = 0; j < nevLocal; ++j)
        memcpy(S+j*numEigen, S+j*localSize, numEigen*sizeof(double));
      localSize = numEigen;
    }

    if ((localSize == numEigen+2*workingSize) && (nevLocal <= numEigen+workingSize)) {
      for (j = 0; j < nevLocal; ++j)
        memcpy(S+j*(numEigen+workingSize), S+j*localSize, (numEigen+workingSize)*sizeof(double));
      localSize = numEigen + workingSize;
    }

    // Update the spaces
    // Note: Use R as a temporary work space
    timeLocalUpdate -= MyWatch.WallTime();

    memcpy(R.Values(), X.Values(), xr*numEigen*sizeof(double));
    callBLAS.GEMM('N', 'N', xr, numEigen, numEigen, 1.0, R.Values(), xr, S, localSize,
                  0.0, X.Values(), xr);
    if (M) {
      memcpy(R.Values(), MX.Values(), xr*numEigen*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, numEigen, numEigen, 1.0, R.Values(), xr, S, localSize,
                    0.0, MX.Values(), xr);
    }
    memcpy(R.Values(), KX.Values(), xr*numEigen*sizeof(double));
    callBLAS.GEMM('N', 'N', xr, numEigen, numEigen, 1.0, R.Values(), xr, S, localSize,
                  0.0, KX.Values(), xr);

    if (localSize == numEigen + workingSize) {
      callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, HI.Values(), xr,
                    S + numEigen, localSize, 0.0, P.Values(), xr);
      if (M) {
        callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, MHI.Values(), xr,
                      S + numEigen, localSize, 0.0, MP.Values(), xr);
      }
      callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, KHI.Values(), xr,
                    S + numEigen, localSize, 0.0, KP.Values(), xr);
    }

    if (localSize > numEigen + workingSize) {
      memcpy(RI.Values(), PI.Values(), xr*workingSize*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, HI.Values(), xr,
                    S + numEigen, localSize, 0.0, P.Values(), xr);
      callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, RI.Values(), xr,
                    S + numEigen + workingSize, localSize, 1.0, P.Values(), xr);
      if (M) {
        memcpy(RI.Values(), MPI.Values(), xr*workingSize*sizeof(double));
        callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, MHI.Values(), xr,
                      S + numEigen, localSize, 0.0, MP.Values(), xr);
        callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, RI.Values(), xr,
                      S + numEigen + workingSize, localSize, 1.0, MP.Values(), xr);
      } 
      memcpy(RI.Values(), KPI.Values(), xr*workingSize*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, KHI.Values(), xr,
                    S + numEigen, localSize, 0.0, KP.Values(), xr);
      callBLAS.GEMM('N', 'N', xr, numEigen, workingSize, 1.0, RI.Values(), xr,
                    S + numEigen + workingSize, localSize, 1.0, KP.Values(), xr);
    }

    if (localSize > numEigen) {
      callBLAS.AXPY(xr*numEigen, 1.0, P.Values(), X.Values());
      if (M)
        callBLAS.AXPY(xr*numEigen, 1.0, MP.Values(), MX.Values());
      callBLAS.AXPY(xr*numEigen, 1.0, KP.Values(), KX.Values());
    }
    timeLocalUpdate += MyWatch.WallTime();

    // Compute the residuals
    timeResidual -= MyWatch.WallTime();
    memcpy(R.Values(), KX.Values(), xr*numEigen*sizeof(double));
    for (j = 0; j < numEigen; ++j) {
      callBLAS.AXPY(xr, -theta[j], MX.Values() + j*xr, R.Values() + j*xr);
    }
    timeResidual += MyWatch.WallTime();
    residual += numEigen;

    // Compute the norms of the residuals
    timeNorm -= MyWatch.WallTime();
    if (vectWeight) {
      R.NormWeighted(*vectWeight, normR);
    }
    else {
      R.Norm2(normR);
    }
    // Scale the norms of residuals with the eigenvalues
    for (j = 0; j < numEigen; ++j) {
      normR[j] = (theta[j] == 0.0) ? normR[j] : normR[j]/theta[j];
    }
    timeNorm += MyWatch.WallTime();

    // When required, monitor some orthogonalities
    if (verbose > 2) {
      accuracyCheck(&X, &MX, &R);
    } // if (verbose > 2)

    if (localSize == numEigen + workingSize)
      localSize += workingSize;

    // Count the converged eigenvectors
    timeNorm -= MyWatch.WallTime();
    knownEV = 0;
    for (i=0; i < numEigen; ++i) {
      if (normR[i] < tolEigenSolve) {
        lambda[i] = theta[i];
        knownEV += 1;
      }
    }
    timeNorm += MyWatch.WallTime();

    // Store the residual history
    if (localVerbose > 2) {
      memcpy(resHistory + historyCount, normR, numEigen*sizeof(double));
      historyCount += numEigen;
    }

    // Print information on current iteration
    if (localVerbose > 0) {
      cout << " Iteration " << outerIter;
      cout << " - Number of converged eigenvectors " << knownEV << endl;
    }

    if (localVerbose > 1) {
      cout << endl;
      cout.precision(2);
      cout.setf(ios::scientific, ios::floatfield);
      for (i=0; i<numEigen; ++i) {
        cout << " Iteration " << outerIter << " - Scaled Norm of Residual " << i;
        cout << " = " << normR[i] << endl;
      }
      cout << endl;
      cout.precision(2);
      for (i=0; i<numEigen; ++i) {
        cout << " Iteration " << outerIter << " - Ritz eigenvalue " << i;
        cout.setf((fabs(theta[i]) < 0.01) ? ios::scientific : ios::fixed, ios::floatfield);
        cout << " = " << theta[i] << endl;
      }
      cout << endl;
    }

    // Convergence test
    if (knownEV >= numEigen) {
      if (localVerbose == 1) {
        cout << endl;
        cout.precision(2);
        cout.setf(ios::scientific, ios::floatfield);
        for (i=0; i<numEigen; ++i) {
          cout << " Iteration " << outerIter << " - Scaled Norm of Residual " << i;
          cout << " = " << normR[i] << endl;
        }
        cout << endl;
      }
      break;
    }

    // Update the sizes
    if ((knownEV > 0) && (localSize > numEigen)) {
      if (localSize == numEigen + workingSize)
        localSize = 2*numEigen - knownEV;
      if (localSize == numEigen + 2*workingSize)
        localSize = 3*numEigen - 2*knownEV;
    }

    // Put the unconverged vectors at the end if needed

    for (i=0; i<numEigen; ++i) {
      if (normR[i] >= tolEigenSolve) {
        firstIndex = i;
        break;
      }
    }

    // Continue the loop if no motion of vectors is necessary
    if (firstIndex == numEigen-1)
      continue;

    while (firstIndex < knownEV) {

      for (j = firstIndex; j < numEigen; ++j) {
        if (normR[j] < tolEigenSolve) {
          callFortran.SWAP(1, normR + j, 1, normR + firstIndex, 1);
          callFortran.SWAP(xr, X.Values() + j*xr, 1, X.Values() + firstIndex*xr, 1);
          callFortran.SWAP(xr, KX.Values() + j*xr, 1, KX.Values() + firstIndex*xr, 1);
          callFortran.SWAP(xr, R.Values() + j*xr, 1, R.Values() + firstIndex*xr, 1);
          callFortran.SWAP(xr, H.Values() + j*xr, 1, H.Values() + firstIndex*xr, 1);
          callFortran.SWAP(xr, KH.Values() + j*xr, 1, KH.Values() + firstIndex*xr, 1);
          callFortran.SWAP(xr, P.Values() + j*xr, 1, P.Values() + firstIndex*xr, 1);
          callFortran.SWAP(xr, KP.Values() + j*xr, 1, KP.Values() + firstIndex*xr, 1);
          if (M) {
            callFortran.SWAP(xr, MX.Values() + j*xr, 1, MX.Values() + firstIndex*xr, 1);
            callFortran.SWAP(xr, MH.Values() + j*xr, 1, MH.Values() + firstIndex*xr, 1);
            callFortran.SWAP(xr, MP.Values() + j*xr, 1, MP.Values() + firstIndex*xr, 1);
          }
          break;
        }
      }

      for (i = firstIndex; i < numEigen; ++i) {
        if (normR[i] >= tolEigenSolve) {
          firstIndex = i;
          break;
        }
      }

    }

  } // for (outerIter = 1; outerIter <= maxIterEigenSolve; ++outerIter)
  timeOuterLoop += MyWatch.WallTime();
  highMem = (highMem > currentSize()) ? highMem : currentSize();

  // Clean memory
  delete[] work1;
  delete[] work2;
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


void KnyazevLOBPCG::accuracyCheck(const Epetra_MultiVector *X, const Epetra_MultiVector *MX,
                       const Epetra_MultiVector *R) const {

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

  if ((R) && (X)) {
    tmp = myVerify.errorOrthogonality(X, R);
    if (myPid == 0)
      cout << " >> Orthogonality X^T R up to " << tmp << endl;
  }

}


void KnyazevLOBPCG::initializeCounters() {

  historyCount = 0;
  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

  memRequested = 0.0;
  highMem = 0.0;

  massOp = 0;
  numRestart = 0;
  outerIter = 0;
  precOp = 0;
  residual = 0;
  stifOp = 0;

  timeLocalProj = 0.0;
  timeLocalSolve = 0.0;
  timeLocalUpdate = 0.0;
  timeMassOp = 0.0;
  timeNorm = 0.0;
  timeOuterLoop = 0.0;
  timePostProce = 0.0;
  timePrecOp = 0.0;
  timeResidual = 0.0;
  timeStifOp = 0.0;


}


void KnyazevLOBPCG::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: LOBPCG (Knyazev's version) with Cholesky-based local eigensolver\n";
    cout << " Block Size: " << blockSize << endl;
  }

}


void KnyazevLOBPCG::historyInfo() const {

  if (resHistory) {
    int j;
    cout << " Block Size: " << blockSize << endl;
    cout << endl;
    cout << " Residuals\n";
    cout << endl;
    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);
    for (j = 0; j < historyCount; ++j) {
      cout << resHistory[j] << endl;
    }
    cout << endl;
  }

}


void KnyazevLOBPCG::memoryInfo() const {

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


void KnyazevLOBPCG::operationInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Operations ---\n\n";
    cout << " Total number of mass matrix multiplications      = ";
    cout.width(9);
    cout << massOp << endl;
    cout << " Total number of stiffness matrix operations      = ";
    cout.width(9);
    cout << stifOp << endl;
    cout << " Total number of preconditioner applications      = ";
    cout.width(9);
    cout << precOp << endl;
    cout << " Total number of computed eigen-residuals         = ";
    cout.width(9);
    cout << residual << endl;
    cout << "\n";
    cout << " Total number of outer iterations                 = ";
    cout.width(9);
    cout << outerIter << endl;
    cout << "       Number of restarts                         = ";
    cout.width(9);
    cout << numRestart << endl;
    cout << "\n";
  }

}


void KnyazevLOBPCG::timeInfo() const {

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
    cout << "       Time for residual computations            = ";
    cout.width(9);
    cout << timeResidual << " s     ";
    cout.width(5);
    cout << 100*timeResidual/timeOuterLoop << " %\n";
    cout << "       Time for residuals norm computations      = ";
    cout.width(9);
    cout << timeNorm << " s     ";
    cout.width(5);
    cout << 100*timeNorm/timeOuterLoop << " %\n";
    cout << "\n";
    cout << " Total time for post-processing                  = ";
    cout.width(9);
    cout << timePostProce << " s\n";
    cout << endl;
  } // if (myPid == 0)

}


