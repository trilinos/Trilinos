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
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#include "Davidson.h"


Davidson::Davidson(const Epetra_Comm &_Comm, const Epetra_Operator *KK, 
                     const Epetra_Operator *PP, int _blk, int _numBlk,
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
           numBlock(_numBlk),
           normWeight(0),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           maxSpaceSize(0),
           sumSpaceSize(0),
           spaceSizeHistory(0),
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


Davidson::Davidson(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     const Epetra_Operator *MM, const Epetra_Operator *PP, int _blk, int _numBlk,
                     double _tol, int _maxIter, int _verb, double *_weight) :
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
           numBlock(_numBlk),
           normWeight(_weight),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           maxSpaceSize(0),
           sumSpaceSize(0),
           spaceSizeHistory(0),
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


Davidson::~Davidson() {

  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

  if (spaceSizeHistory) {
    delete[] spaceSizeHistory;
    spaceSizeHistory = 0;
  }

}


int Davidson::solve(int numEigen, Epetra_MultiVector &Q, double *lambda) {

  return Davidson::reSolve(numEigen, Q, lambda);

}


int Davidson::reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV) {

  // Computes the smallest eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = M X Lambda
  // 
  // using a generalized Davidson algorithm
  //
  // Note that if M is not specified, then  K X = X Lambda is solved.
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
  // MX = Working vectors (storing M*X if M is specified, else pointing to X)
  // KX = Working vectors (storing K*X)
  //
  // R = Residuals

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

  int tmp;
  tmp = (M == 0) ? 2*blockSize*xr : 3*blockSize*xr;

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

  Epetra_MultiVector KX(View, Q.Map(), tmpD, xr, blockSize);
  tmpD = tmpD + xr*blockSize;

  Epetra_MultiVector MX(View, Q.Map(), (M) ? tmpD : X.Values(), xr, blockSize);
  tmpD = (M) ? tmpD + xr*blockSize : tmpD;

  Epetra_MultiVector R(View, Q.Map(), tmpD, xr, blockSize);

  // Define arrays
  //
  // theta = Store the local eigenvalues (size: dimSearch)
  // normR = Store the norm of residuals (size: blockSize)
  //
  // KK = Local stiffness matrix         (size: dimSearch x dimSearch)
  //
  // S = Local eigenvectors              (size: dimSearch x dimSearch)
  //
  // tmpKK = Local workspace             (size: blockSize x blockSize)

  int lwork2 = blockSize + dimSearch + 2*dimSearch*dimSearch + blockSize*blockSize;
  double *work2 = new (nothrow) double[lwork2];
  if (work2 == 0) {
    if (vectWeight)
      delete vectWeight;
    delete[] work1;
    info = -30;
    return info;
  }

  memRequested += sizeof(double)*lwork2/(1024.0*1024.0);
  highMem = (highMem > currentSize()) ? highMem : currentSize();

  tmpD = work2;

  double *theta = tmpD;
  tmpD = tmpD + dimSearch;

  double *normR = tmpD;
  tmpD = tmpD + blockSize;

  double *KK = tmpD;
  tmpD = tmpD + dimSearch*dimSearch;
  memset(KK, 0, dimSearch*dimSearch*sizeof(double));

  double *S = tmpD;
  tmpD = tmpD + dimSearch*dimSearch;

  double *tmpKK = tmpD;

  // Define an array to store the residuals history
  if (localVerbose > 2) {
    resHistory = new (nothrow) double[maxIterEigenSolve*blockSize];
    spaceSizeHistory = new (nothrow) int[maxIterEigenSolve];
    if ((resHistory == 0) || (spaceSizeHistory == 0)) {
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

  bool criticalExit = false;

  int bStart = 0;
  int offSet = 0;
  numBlock = (dimSearch/blockSize) - (knownEV/blockSize);

  int nFound = blockSize;
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
    cout << " *|* Algorithm = Davidson algorithm (block version)" << endl;
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
            Epetra_MultiVector copyQ(View, X, 0, knownEV + localSize);
            info = modalTool.massOrthonormalize(Xcurrent, MX, M, copyQ, nFound, 0, R.Values());
          }
        }
        nFound = 0;
      }
      else {
        Epetra_MultiVector copyQ(View, X, 0, knownEV + localSize);
        info = modalTool.massOrthonormalize(Xcurrent, MX, M, copyQ, blockSize, 0, R.Values());
      }
      timeOrtho += MyWatch.WallTime();

      // Exit the code when the number of vectors exceeds the space dimension
      if (info < 0) {
        delete[] work1;
        delete[] work2;
        if (vectWeight)
          delete vectWeight;
        return -10;
      }

      timeStifOp -= MyWatch.WallTime();
      K->Apply(Xcurrent, KX);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;

      // Check the orthogonality properties of X
      if (verbose > 2) {
        if (knownEV + localSize == 0)
          accuracyCheck(&Xcurrent, &MX, 0);
        else {
          Epetra_MultiVector copyQ(View, X, 0, knownEV + localSize);
          accuracyCheck(&Xcurrent, &MX, &copyQ);
        }
        if (localVerbose > 0)
          cout << endl;
      } // if (verbose > 2)

      // Define the local stiffness matrix
      // Note: S is used as a workspace
      timeLocalProj -= MyWatch.WallTime();
      for (j = 0; j <= nb; ++j) {
        callBLAS.GEMM('T', 'N', blockSize, blockSize, xr,
                      1.0, X.Values()+(knownEV+j*blockSize)*xr, xr, KX.Values(), xr,
                      0.0, tmpKK, blockSize);
        MyComm.SumAll(tmpKK, S, blockSize*blockSize);
        int iC;
        for (iC = 0; iC < blockSize; ++iC) {
          double *Kpointer = KK + localSize*dimSearch + j*blockSize + iC*dimSearch;
          memcpy(Kpointer, S + iC*blockSize, blockSize*sizeof(double));
        }
      }
      timeLocalProj += MyWatch.WallTime();

      // Perform a spectral decomposition
      timeLocalSolve -= MyWatch.WallTime();
      int nevLocal = localSize + blockSize;
      info = modalTool.directSolver(localSize+blockSize, KK, dimSearch, 0, 0,
                                    nevLocal, S, dimSearch, theta, localVerbose, 10);
      timeLocalSolve += MyWatch.WallTime();

      if (info != 0) {
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

      // Update the search space
      // Note: Use KX as a workspace
      timeLocalUpdate -= MyWatch.WallTime();
      callBLAS.GEMM('N', 'N', xr, blockSize, localSize+blockSize, 1.0, X.Values()+knownEV*xr, xr,
                    S, dimSearch, 0.0, KX.Values(), xr);
      timeLocalUpdate += MyWatch.WallTime();

      // Apply the mass matrix for the next block
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(KX, MX);
      timeMassOp += MyWatch.WallTime();
      massOp += blockSize;

      // Apply the stiffness matrix for the next block
      timeStifOp -= MyWatch.WallTime();
      K->Apply(KX, R);
      timeStifOp += MyWatch.WallTime();
      stifOp += blockSize;

      // Form the residuals
      timeResidual -= MyWatch.WallTime();
      if (M) {
        for (j = 0; j < blockSize; ++j) {
          callBLAS.AXPY(xr, -theta[j], MX.Values() + j*xr, R.Values() + j*xr);
        }
      }
      else {
        // Note KX contains the updated block
        for (j = 0; j < blockSize; ++j) {
          callBLAS.AXPY(xr, -theta[j], KX.Values() + j*xr, R.Values() + j*xr);
        }
      }
      timeResidual += MyWatch.WallTime();
      residual += blockSize;

      // Compute the norm of residuals
      timeNorm -= MyWatch.WallTime();
      if (vectWeight) {
        R.NormWeighted(*vectWeight, normR);
      }
      else {
        R.Norm2(normR);
      }
      // Scale the norms of residuals with the eigenvalues
      // Count the number of converged eigenvectors
      nFound = 0;
      for (j = 0; j < blockSize; ++j) {
        normR[j] = (theta[j] == 0.0) ? normR[j] : normR[j]/theta[j];
        if (normR[j] < tolEigenSolve)
          nFound += 1;
      } // for (j = 0; j < blockSize; ++j)
      timeNorm += MyWatch.WallTime();

      // Store the residual history
      if (localVerbose > 2) {
        memcpy(resHistory + historyCount*blockSize, normR, blockSize*sizeof(double));
        spaceSizeHistory[historyCount] = localSize + blockSize;
        historyCount += 1;
      }
      maxSpaceSize = (maxSpaceSize > localSize+blockSize) ? maxSpaceSize : localSize+blockSize;
      sumSpaceSize += localSize + blockSize;

      // Print information on current iteration
      if (localVerbose > 0) {
        cout << " Iteration " << outerIter << " - Number of converged eigenvectors ";
        cout << knownEV + nFound << endl;
      } // if (localVerbose > 0)

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
        for (i=0; i<nevLocal; ++i) {
          cout << " Iteration " << outerIter << " - Ritz eigenvalue " << i;
          cout.setf((fabs(theta[i]) < 0.01) ? ios::scientific : ios::fixed, ios::floatfield);  
          cout << " = " << theta[i] << endl;
        }
        cout << endl;
      }

      // Exit the loop to treat the converged eigenvectors
      if (nFound > 0) {
        nb += 1;
        offSet = 0;
        break;
      }

      // Apply the preconditioner on the residuals
      // Note: Use KX as a workspace
      if (maxBlock == 1) {
        if (Prec) {
          timePrecOp -= MyWatch.WallTime();
          Prec->ApplyInverse(R, Xcurrent);
          timePrecOp += MyWatch.WallTime();
          precOp += blockSize;
        }
        else {
          memcpy(Xcurrent.Values(), R.Values(), blockSize*xr*sizeof(double));
        }
        timeRestart -= MyWatch.WallTime();
        Xcurrent.Update(1.0, KX, -1.0);
        timeRestart += MyWatch.WallTime();
        break;
      } // if (maxBlock == 1)

      if (nb == maxBlock - 1) {
        nb += 1;
        break;
      }

      Epetra_MultiVector Xnext(View, X, knownEV+localSize+blockSize, blockSize);
      if (Prec) {
        timePrecOp -= MyWatch.WallTime();
        Prec->ApplyInverse(R, Xnext);
        timePrecOp += MyWatch.WallTime();
        precOp += blockSize;
      }
      else {
        memcpy(Xnext.Values(), R.Values(), blockSize*xr*sizeof(double));
      }

    } // for (nb = bStart; nb < maxBlock; ++nb)

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

    // Treat the particular case of 1 block
    if (maxBlock == 1) {
      if (nFound > 0) {
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
        Epetra_MultiVector Xnext(View, X, knownEV + blockSize - nFound, nFound);
        Xnext.Random();
      }
      else {
        nFound = blockSize;
      }
      continue;
    }

    // Define the restarting block when maxBlock > 1
    if (nFound > 0) {
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

    } // if (nFound > 0)

    // Define the restarting size
    bStart = ((nb - offSet) > 2) ? (nb - offSet)/2 : 0;

    // Define the restarting space and local stiffness
    timeRestart -= MyWatch.WallTime();
    memset(KK, 0, nb*blockSize*dimSearch*sizeof(double));
    for (j = 0; j < bStart*blockSize; ++j) {
      KK[j + j*dimSearch] = theta[j + nFound];
    }
    // Form the restarting space
    int oldCol = nb*blockSize;
    int newCol = nFound + (bStart+1)*blockSize;
    newCol = (newCol > oldCol) ? oldCol : newCol;
    callFortran.GEQRF(oldCol, newCol, S, dimSearch, theta, R.Values(), xr*blockSize, &info);
    callFortran.ORMQR('R', 'N', xr, oldCol, newCol, S, dimSearch, theta,
                      X.Values()+knownEV*xr, xr, R.Values(), blockSize*xr, &info);
    timeRestart += MyWatch.WallTime();

    if (nFound == 0)
      offSet += 1;

    knownEV += nFound;
    maxBlock = (dimSearch/blockSize) - (knownEV/blockSize);

    // Put random vectors if the Rayleigh Ritz vectors are not enough
    newCol = nFound + (bStart+1)*blockSize;
    if (newCol > oldCol) {
      Epetra_MultiVector Xnext(View, X, knownEV+blockSize-nFound, nFound);
      Xnext.Random();
      continue;
    }

    nFound = 0;

  } // while (outerIter <= maxIterEigenSolve)
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


void Davidson::accuracyCheck(const Epetra_MultiVector *X, const Epetra_MultiVector *MX,
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


int Davidson::minimumSpaceDimension(int nev) const {

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


void Davidson::initializeCounters() {

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


void Davidson::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: Davidson algorithm (block version)\n";
    cout << " Block Size: " << blockSize << endl;
    cout << " Number of Blocks kept: " << numBlock << endl;
  }

}


void Davidson::historyInfo() const {

  if (resHistory) {
    int j;
    cout << " Block Size: " << blockSize << endl;
    cout << endl;
    cout << " Residuals    Search Space Dim.\n";
    cout << endl;
    cout.precision(4);
    cout.setf(ios::scientific, ios::floatfield);
    for (j = 0; j < historyCount; ++j) {
      int ii;
      for (ii = 0; ii < blockSize; ++ii) {
        cout << resHistory[blockSize*j + ii] << "      ";
        cout.width(4);
        cout << spaceSizeHistory[j] << endl;
      }
    }
    cout << endl;
  }

}


void Davidson::memoryInfo() const {

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


void Davidson::operationInfo() const {

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


void Davidson::timeInfo() const {

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
//    cout << "                 Q^T mult.  :: ";
//    cout.width(9);
//    cout << modalTool.getTimeProj_QtMult() << " s\n";
//    cout << "                 Q mult.    :: ";
//    cout.width(9);
//    cout << modalTool.getTimeProj_QMult() << " s\n";
//    cout << "                 Mass mult. :: ";
//    cout.width(9);
//    cout << modalTool.getTimeProj_MassMult() << " s\n";
    cout << "            Normalization step       : ";
    cout.width(9);
    cout << modalTool.getTimeNorm() << " s\n";
//    cout << "                 Mass mult. :: ";
//    cout.width(9);
//    cout << modalTool.getTimeNorm_MassMult() << " s\n";
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


