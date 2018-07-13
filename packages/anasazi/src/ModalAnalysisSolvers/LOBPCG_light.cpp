// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
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

#include "LOBPCG_light.h"


LOBPCG_light::LOBPCG_light(const Epetra_Comm &_Comm, const Epetra_Operator *KK, 
                     const Epetra_Operator *PP, int _blk,
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
           blockSize(_blk),
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
           timeOrtho(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timePrecOp(0.0),
           timeResidual(0.0),
           timeRestart(0.0),
           timeStifOp(0.0)
           {

}


LOBPCG_light::LOBPCG_light(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     const Epetra_Operator *MM, const Epetra_Operator *PP, int _blk,
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
           blockSize(_blk),
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
           timeOrtho(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timePrecOp(0.0),
           timeResidual(0.0),
           timeRestart(0.0),
           timeStifOp(0.0)
           {

}


LOBPCG_light::~LOBPCG_light() {

  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

}


int LOBPCG_light::reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV) {

  // Computes the smallest eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = M X Lambda
  // 
  // using a modified Locally Optimal Block Preconditioned Conjugate
  // Gradient method. 
  //
  // Note that if M is not specified, then  K X = X Lambda is solved.
  // 
  // Ref: A. Knyazev, "Toward the optimal preconditioned eigensolver:
  // Locally optimal block preconditioner conjugate gradient method",
  // SIAM J. Sci. Comput., vol 23, n 2, pp. 517-541
  // 
  // Input variables:
  // 
  // numEigen  (integer) = Number of eigenmodes requested
  // 
  // Q (Epetra_MultiVector) = Converged eigenvectors
  //                   The number of columns of Q must be equal to numEigen + blockSize.
  //                   The rows of Q are distributed across processors.
  //                   At exit, the first numEigen columns contain the eigenvectors requested.
  // 
  // lambda (array of doubles) = Converged eigenvalues
  //                   At input, it must be of size numEigen + blockSize.
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
    return numEigen;
  }

  int info = myVerify.inputArguments(numEigen, K, M, Prec, Q, minimumSpaceDimension(numEigen));
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
  // R = Residuals
  //
  // H = Preconditioned search space
  //
  // P = Search directions

  int xr = Q.MyLength();
  Epetra_MultiVector X(View, Q, numEigen, blockSize);
  X.Random();

  int tmp;
  tmp = 3*blockSize*xr;

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

  Epetra_MultiVector R(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + xr*blockSize;

  Epetra_MultiVector H(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + xr*blockSize;

  Epetra_MultiVector P(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + xr*blockSize;

  // Define arrays
  //
  // theta = Store the local eigenvalues (size: 3*blockSize)
  // normR = Store the norm of residuals (size: blockSize)
  //
  // MM = Local mass matrix              (size: 3*blockSize x 3*blockSize)
  // KK = Local stiffness matrix         (size: 3*blockSize x 3*blockSize)
  //
  // S = Local eigenvectors              (size: 3*blockSize x 3*blockSize)

  int lwork2;
  lwork2 = 4*blockSize + 27*blockSize*blockSize;
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
  tmpD = tmpD + 3*blockSize;

  double *normR = tmpD;
  tmpD = tmpD + blockSize;

  double *MM = tmpD;
  tmpD = tmpD + 9*blockSize*blockSize;

  double *KK = tmpD;
  tmpD = tmpD + 9*blockSize*blockSize;

  double *S = tmpD;

  memRequested += sizeof(double)*lwork2/(1024.0*1024.0);

  // Define an array to store the residuals history
  if (localVerbose > 2) {
    resHistory = new (nothrow) double[maxIterEigenSolve*blockSize];
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
  int twoBlocks = 2*blockSize;
  int threeBlocks = 3*blockSize;
  int nFound = blockSize;
  int i, j;

  double minNorm = 1000*tolEigenSolve;
  bool hasNextVectors = false;

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
    cout << " *|* Algorithm = LOBPCG (Small memory requirements)" << endl;
    cout << " *|* Size of blocks = " << blockSize << endl;
    cout << " *|* Number of requested eigenvalues = " << numEigen << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " *|* Tolerance for convergence = " << tolEigenSolve << endl;
    cout << " *|* Norm used for convergence: ";
    if (normWeight){
      cout << "weighted L2-norm with user-provided weights" << endl;
    }
    else {
      cout << "L^2-norm" << endl;
    }
    if (startingEV > 0)
      cout << " *|* Input converged eigenvectors = " << startingEV << endl;
    cout << "\n -- Start iterations -- \n";
  }

  timeOuterLoop -= MyWatch.WallTime();
  for (outerIter = 1; outerIter <= maxIterEigenSolve; ++outerIter) {

    highMem = (highMem > currentSize()) ? highMem : currentSize();

    if ((outerIter == 1) || (reStart == true)) {

      reStart = false;
      localSize = blockSize;

      hasNextVectors = false;
      minNorm = 1000*tolEigenSolve;

      // Apply the mass matrix to X
      // Note: Use P as work space
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(X, P);
      else
        memcpy(P.Values(), X.Values(), xr*blockSize*sizeof(double));
      timeMassOp += MyWatch.WallTime();
      massOp += blockSize;

      if ((knownEV > 0) && (nFound > 0)) {
        // Orthonormalize X against the known eigenvectors with Gram-Schmidt
        // Note: Use R as a temporary work space
        Epetra_MultiVector copyQ(View, Q, 0, knownEV);
        timeOrtho -= MyWatch.WallTime();
        info = modalTool.massOrthonormalize(X, P, M, copyQ, nFound, 1, R.Values());
        timeOrtho += MyWatch.WallTime();
        // Exit the code if the orthogonalization did not succeed
        if (info < 0) {
          info = -10;
          if (vectWeight)
            delete vectWeight;
          delete[] work1;
          delete[] work2;
          return info;
        }
      } // if ((knownEV > 0) && (nFound > 0))

    } // if ((outerIter == 1) || (reStart == true))
    else {

      // Apply the preconditioner on the residuals
      if (Prec) {
        timePrecOp -= MyWatch.WallTime();
        Prec->ApplyInverse(R, H);
        timePrecOp += MyWatch.WallTime();
        precOp += blockSize;
      }
      else {
        memcpy(H.Values(), R.Values(), xr*blockSize*sizeof(double));
      }

      // Apply the mass matrix on H
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(H, R);
      else
        memcpy(R.Values(), H.Values(), xr*blockSize*sizeof(double));
      timeMassOp += MyWatch.WallTime();
      massOp += blockSize;

      if (knownEV > 0) {
        // Orthogonalize H against the known eigenvectors
        Epetra_MultiVector copyQ(View, Q, 0, knownEV);
        timeOrtho -= MyWatch.WallTime();
        modalTool.massOrthonormalize(H, R, M, copyQ, blockSize, 1);
        timeOrtho += MyWatch.WallTime();
      }

      if (localSize == blockSize)
        localSize += blockSize;

    } // if ((outerIter == 1) || (reStart==true))

    // Form "local" mass and stiffness matrices
    // Note: Use S as a temporary workspace
    if (localSize == blockSize) {
      // P stores M*X
      timeLocalProj -= MyWatch.WallTime();
      modalTool.localProjection(blockSize, blockSize, xr, X.Values(), xr, P.Values(), xr,
                                MM, localSize, S);
      timeLocalProj += MyWatch.WallTime();
      // Put K*X in H
      timeStifOp -= MyWatch.WallTime();
      K->Apply(X, H);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;
      timeLocalProj -= MyWatch.WallTime();
      modalTool.localProjection(blockSize, blockSize, xr, X.Values(), xr, H.Values(), xr,
                                KK, localSize, S);
      timeLocalProj += MyWatch.WallTime();
    }
    else {
      // Put diagonal matrices in first block
      timeLocalProj -= MyWatch.WallTime();
      for (i = 0; i < blockSize; ++i) {
        memset(KK + i*localSize, 0, i*sizeof(double));
        KK[i + i*localSize] = theta[i];
        memset(MM + i*localSize, 0, i*sizeof(double));
        MM[i + i*localSize] = 1.0;
      }
      timeLocalProj += MyWatch.WallTime();
      // R stores M*H
      timeLocalProj -= MyWatch.WallTime();
      modalTool.localProjection(blockSize, blockSize, xr, X.Values(), xr, R.Values(), xr,
                                MM + blockSize*localSize, localSize, S);
      modalTool.localProjection(blockSize, blockSize, xr, H.Values(), xr, R.Values(), xr,
                                MM + blockSize*localSize + blockSize, localSize, S);
      timeLocalProj += MyWatch.WallTime();
      // Put K*H in R
      timeStifOp -= MyWatch.WallTime();
      K->Apply(H, R);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;
      timeLocalProj -= MyWatch.WallTime();
      modalTool.localProjection(blockSize, blockSize, xr, X.Values(), xr, R.Values(), xr,
                                KK + blockSize*localSize, localSize, S);
      modalTool.localProjection(blockSize, blockSize, xr, H.Values(), xr, R.Values(), xr,
                                KK + blockSize*localSize + blockSize, localSize, S);
      timeLocalProj += MyWatch.WallTime();
      if (localSize > twoBlocks) {
        // Put M*P in R
        timeMassOp -= MyWatch.WallTime();
        if (M)
          M->Apply(P, R);
        else
          memcpy(R.Values(), P.Values(), xr*blockSize*sizeof(double));
        timeMassOp += MyWatch.WallTime();
        massOp += blockSize;
        timeLocalProj -= MyWatch.WallTime();
        modalTool.localProjection(blockSize, blockSize, xr, X.Values(), xr, R.Values(), xr,
                                  MM + twoBlocks*localSize, localSize, S);
        modalTool.localProjection(blockSize, blockSize, xr, H.Values(), xr, R.Values(), xr,
                                  MM + twoBlocks*localSize + blockSize, localSize, S);
        modalTool.localProjection(blockSize, blockSize, xr, P.Values(), xr, R.Values(), xr,
                                  MM + twoBlocks*localSize + twoBlocks, localSize, S);
        timeLocalProj += MyWatch.WallTime();
        // Put K*P in R
        timeStifOp -= MyWatch.WallTime();
        K->Apply(P, R);
        timeStifOp += MyWatch.WallTime();
        stifOp += blockSize;
        timeLocalProj -= MyWatch.WallTime();
        modalTool.localProjection(blockSize, blockSize, xr, X.Values(), xr, R.Values(), xr,
                                  KK + twoBlocks*localSize, localSize, S);
        modalTool.localProjection(blockSize, blockSize, xr, H.Values(), xr, R.Values(), xr,
                                  KK + twoBlocks*localSize + blockSize, localSize, S);
        modalTool.localProjection(blockSize, blockSize, xr, P.Values(), xr, R.Values(), xr,
                                  KK + twoBlocks*localSize + twoBlocks, localSize, S);
        timeLocalProj += MyWatch.WallTime();
      } // if (localSize > twoBlocks)
    } // if (localSize == blockSize)

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
    if ((theta[0] < 0.0) || (nevLocal < blockSize)) {
      if (localVerbose > 0) {
        cout << " Iteration " << outerIter;
        cout << "- Failure for spectral decomposition - RESTART with new random search\n";
      }
      if (blockSize == 1) {
        X.Random();
        nFound = blockSize;
      }
      else {
        Epetra_MultiVector Xinit(View, X, 1, blockSize-1);
        Xinit.Random();
        nFound = blockSize - 1;
      } // if (blockSize == 1)
      reStart = true;
      numRestart += 1;
      info = 0;
      continue;
    } // if ((theta[0] < 0.0) || (nevLocal < blockSize))

    if ((localSize == twoBlocks) && (nevLocal == blockSize)) {
      for (j = 0; j < nevLocal; ++j) 
        memcpy(S + j*blockSize, S + j*twoBlocks, blockSize*sizeof(double)); 
      localSize = blockSize;
    }

    if ((localSize == threeBlocks) && (nevLocal <= twoBlocks)) {
      for (j = 0; j < nevLocal; ++j) 
        memcpy(S + j*twoBlocks, S + j*threeBlocks, twoBlocks*sizeof(double)); 
      localSize = twoBlocks;
    }

    // Update the iterate and, if needed, define the next Ritz vectors
    if (localSize == twoBlocks) {
      if (minNorm < 10*tolEigenSolve) {
        timeRestart -= MyWatch.WallTime();
        int numVec = (numEigen > knownEV + blockSize) ? blockSize : numEigen - knownEV - 1;
        double *pointerS = S + blockSize*localSize;
        double *pointerQ = Q.Values() + knownEV*xr;
        callBLAS.GEMM('N', 'N', xr, numVec, blockSize, 1.0, X.Values(), xr,
                      pointerS, localSize, 0.0, pointerQ, xr);
        callBLAS.GEMM('N', 'N', xr, numVec, blockSize, 1.0, H.Values(), xr,
                      pointerS + blockSize, localSize, 1.0, pointerQ, xr);
        hasNextVectors = true;
        timeRestart += MyWatch.WallTime();
      }
      timeLocalUpdate -= MyWatch.WallTime();
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, H.Values(), xr,
                    S + blockSize, localSize, 0.0, P.Values(), xr);
      memcpy(R.Values(), X.Values(), xr*blockSize*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, R.Values(), xr,
                    S, localSize, 0.0, X.Values(), xr);
      X.Update(1.0, P, 1.0);
      timeLocalUpdate += MyWatch.WallTime();
    }
    if (localSize == threeBlocks) {
      if (minNorm < 10*tolEigenSolve) {
        timeRestart -= MyWatch.WallTime();
        int numVec = (numEigen > knownEV + blockSize) ? blockSize : numEigen - knownEV - 1;
        double *pointerS = S + blockSize*localSize;
        double *pointerQ = Q.Values() + knownEV*xr;
        callBLAS.GEMM('N', 'N', xr, numVec, blockSize, 1.0, X.Values(), xr,
                      pointerS, localSize, 0.0, pointerQ, xr);
        // Note: We use the contiguity of [H, P] in the memory
        callBLAS.GEMM('N', 'N', xr, numVec, twoBlocks, 1.0, H.Values(), xr,
                      pointerS + blockSize, localSize, 1.0, pointerQ, xr);
        hasNextVectors = true;
        timeRestart += MyWatch.WallTime();
      }
      timeLocalUpdate -= MyWatch.WallTime();
      memcpy(R.Values(), P.Values(), xr*blockSize*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, H.Values(), xr,
                    S + blockSize, localSize, 0.0, P.Values(), xr);
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, R.Values(), xr,
                    S + twoBlocks, localSize, 1.0, P.Values(), xr);
      memcpy(R.Values(), X.Values(), xr*blockSize*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, R.Values(), xr,
                    S, localSize, 0.0, X.Values(), xr);
      X.Update(1.0, P, 1.0);
      timeLocalUpdate += MyWatch.WallTime();
    }

    // Compute the residuals
    if (localSize == blockSize) {
      timeResidual -= MyWatch.WallTime();
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, H.Values(), xr, S, blockSize,
                    0.0, R.Values(), xr);
      timeResidual += MyWatch.WallTime();
      timeLocalUpdate -= MyWatch.WallTime();
      memcpy(H.Values(), X.Values(), xr*blockSize*sizeof(double));
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, 1.0, H.Values(), xr, S, blockSize,
                    0.0, X.Values(), xr);
      timeLocalUpdate += MyWatch.WallTime();
      // Note: We scale S
      timeResidual -= MyWatch.WallTime();
      for (j = 0; j < blockSize; ++j)
        callBLAS.SCAL(localSize, theta[j], S + j*localSize);
      callBLAS.GEMM('N', 'N', xr, blockSize, blockSize, -1.0, P.Values(), xr, S, blockSize,
                    1.0, R.Values(), xr);
      timeResidual += MyWatch.WallTime();
      // Note: S is not scaled back to its original value
    }
    else {
      timeStifOp -= MyWatch.WallTime();
      K->Apply(X, R);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(X, H);
      else
        memcpy(H.Values(), X.Values(), xr*blockSize*sizeof(double));
      timeMassOp += MyWatch.WallTime();
      massOp += blockSize;
      timeResidual -= MyWatch.WallTime();
      for (j = 0; j < blockSize; ++j)
        callBLAS.AXPY(xr, -theta[j], H.Values() + j*xr, R.Values() + j*xr);
      timeResidual += MyWatch.WallTime();
    }

    // Compute the norms of the residuals
    timeNorm -= MyWatch.WallTime();
    if (vectWeight)
      R.NormWeighted(*vectWeight, normR);
    else
      R.Norm2(normR);

    // Scale the norms of residuals with the eigenvalues
    // Count the converged eigenvectors
    nFound = 0;
    for (j = 0; j < blockSize; ++j) {
      normR[j] = (theta[j] == 0.0) ? normR[j] : normR[j]/theta[j];
      if (normR[j] < tolEigenSolve) 
        nFound += 1;
      minNorm = (normR[j] < minNorm) ? normR[j] : minNorm;
    }
    timeNorm += MyWatch.WallTime();

    // Store the residual history
    if (localVerbose > 2) {
      memcpy(resHistory + historyCount*blockSize, normR, blockSize*sizeof(double));
      historyCount += 1;
    }

    // Print information on current iteration
    if (localVerbose > 0) {
      cout << " Iteration " << outerIter << " - Number of converged eigenvectors ";
      cout << knownEV + nFound << endl;
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
      for (i=0; i<blockSize; ++i) {
        cout << " Iteration " << outerIter << " - Ritz eigenvalue " << i;
        cout.setf((fabs(theta[i]) < 0.01) ? ios::scientific : ios::fixed, ios::floatfield);
        cout << " = " << theta[i] << endl;
      }
      cout << endl;
    }

    if (nFound == 0) {
      // When required, monitor some orthogonalities
      if (verbose > 2) {
        if (knownEV == 0) {
          accuracyCheck(&X, 0, &R, 0, 0, 0);
        }
        else {
          Epetra_MultiVector copyQ(View, Q, 0, knownEV);
          accuracyCheck(&X, 0, &R, &copyQ, (localSize>blockSize) ? &H : 0,
                        (localSize>twoBlocks) ? &P : 0);
        }
      } // if (verbose > 2)
      if (localSize < threeBlocks)
        localSize += blockSize;
      continue;
    } // if (nFound == 0)

    // Exit the loop if enough vectors have converged
    if (knownEV + nFound >= numEigen) {
      for (j = 0; j < blockSize; ++j) {
        if (normR[j] < tolEigenSolve) {
          lambda[knownEV] = theta[j];
          memcpy(Q.Values() + knownEV*xr, X.Values() + j*xr, xr*sizeof(double)); 
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
    } // if (knownEV >= numEigen)

    // Store the converged eigenpairs and define the new X iterate
    if (hasNextVectors == true) {
      // Q stores the next Ritz vectors
      for (j = 0; j < blockSize; ++j) {
        if (normR[j] < tolEigenSolve) {
          lambda[knownEV] = theta[j];
          callFortran.SWAP(xr, X.Values() + j*xr, 1, Q.Values() + knownEV*xr, 1);
          knownEV += 1;
        }
      }
      nFound = 0;
    }
    else {
      for (j = 0; j < blockSize; ++j) {
        if (normR[j] < tolEigenSolve) {
          lambda[knownEV] = theta[j];
          memcpy(Q.Values() + knownEV*xr, X.Values() + j*xr, xr*sizeof(double)); 
          knownEV += 1;
        }
      }
      // Replace converged vectors with random vectors
      timeRestart -= MyWatch.WallTime();
      nFound = 0;
      for  (j = 0; j < blockSize; ++j) {
        if (normR[j] >= tolEigenSolve) {
          if (nFound > 0) 
            memcpy(X.Values() + (j-nFound)*xr, X.Values() + j*xr, xr*sizeof(double));
        }
        else 
          nFound += 1;
      }
      if (nFound > 0) {
        // Put new random vectors at the end of the block
        Epetra_MultiVector Xtmp(View, X, blockSize - nFound,  nFound);
        Xtmp.Random();
      }
      timeRestart += MyWatch.WallTime();
    }
    reStart = true;

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


void LOBPCG_light::accuracyCheck(const Epetra_MultiVector *X, const Epetra_MultiVector *MX,
                       const Epetra_MultiVector *R, const Epetra_MultiVector *Q,
                       const Epetra_MultiVector *H, const Epetra_MultiVector *P) const {

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
    if (H) {
      tmp = myVerify.errorOrthogonality(Q, H, M);
      if (myPid == 0)
        cout << " >> Orthogonality Q^T M H up to " << tmp << endl;
    }
    if (P) {
      tmp = myVerify.errorOrthogonality(Q, P, M);
      if (myPid == 0)
        cout << " >> Orthogonality Q^T M P up to " << tmp << endl;
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
    if (H) {
      tmp = myVerify.errorOrthogonality(Q, H, 0);
      if (myPid == 0)
        cout << " >> Orthogonality Q^T H up to " << tmp << endl;
    }
    if (P) {
      tmp = myVerify.errorOrthogonality(Q, P, 0);
      if (myPid == 0)
        cout << " >> Orthogonality Q^T P up to " << tmp << endl;
    }
  }

}


void LOBPCG_light::initializeCounters() {

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
  timeOrtho = 0.0;
  timeOuterLoop = 0.0;
  timePostProce = 0.0;
  timePrecOp = 0.0;
  timeResidual = 0.0;
  timeRestart = 0.0;
  timeStifOp = 0.0;


}


void LOBPCG_light::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: LOBPCG (block version) with storage for 4 blocks of vectors\n";
    cout << " Block Size: " << blockSize << endl;
  }

}


void LOBPCG_light::historyInfo() const {

  if (resHistory) {
    int j;
    cout << " Block Size: " << blockSize << endl;
    cout << endl;
    cout << " Residuals\n";
    cout << endl;
    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);
    for (j = 0; j < historyCount; ++j) {
      int ii;
      for (ii = 0; ii < blockSize; ++ii)
        cout << resHistory[blockSize*j + ii] << endl;
    }
    cout << endl;
  }

}


void LOBPCG_light::memoryInfo() const {

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


void LOBPCG_light::operationInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Operations ---\n\n";
    cout << " Total number of mass matrix multiplications      = ";
    cout.width(9);
    cout << massOp + modalTool.getNumProj_MassMult() + modalTool.getNumNorm_MassMult() << endl;
    cout << "       Mass multiplications in solver loop        = ";
    cout.width(9);
    cout << massOp << endl;
    cout << "       Mass multiplications for orthogonalization = ";
    cout.width(9);
    cout << modalTool.getNumProj_MassMult() << endl;
    cout << "       Mass multiplications for normalization     = ";
    cout.width(9);
    cout << modalTool.getNumNorm_MassMult() << endl;
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


void LOBPCG_light::timeInfo() const {

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
    cout << "                 Q^T mult.  :: ";
    cout.width(9);
    cout << modalTool.getTimeProj_QtMult() << " s\n";
    cout << "                 Q mult.    :: ";
    cout.width(9);
    cout << modalTool.getTimeProj_QMult() << " s\n";
    cout << "                 Mass mult. :: ";
    cout.width(9);
    cout << modalTool.getTimeProj_MassMult() << " s\n";
    cout << "            Normalization step       : ";
    cout.width(9);
    cout << modalTool.getTimeNorm() << " s\n";
    cout << "                 Mass mult. :: ";
    cout.width(9);
    cout << modalTool.getTimeNorm_MassMult() << " s\n";
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
    cout << "       Time for restarting space definition      = ";
    cout.width(9);
    cout << timeRestart << " s     ";
    cout.width(5);
    cout << 100*timeRestart/timeOuterLoop << " %\n";
    cout << "\n";
    cout << " Total time for post-processing                  = ";
    cout.width(9);
    cout << timePostProce << " s\n";
    cout << endl;
  } // if (myPid == 0)

}


