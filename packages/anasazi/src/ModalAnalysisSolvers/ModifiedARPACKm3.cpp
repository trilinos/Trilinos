//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical report SAND2003-1028J.
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

#include "ModifiedARPACKm3.h"


ModifiedARPACKm3::ModifiedARPACKm3(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     double _tol, int _maxIter, int _verb) :
           myVerify(_Comm),
           callBLAS(),
           callFortran(),
           modalTool(_Comm),
           mySort(),
           MyComm(_Comm),
           K(KK),
           M(0),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           normWeight(0),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           numResidual(0),
           orthoOp(0),
           outerIter(0),
           stifOp(0),
           timeMassOp(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timeResidual(0.0),
           timeStifOp(0.0)
           {

}


ModifiedARPACKm3::ModifiedARPACKm3(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     const Epetra_Operator *MM,
                     double _tol, int _maxIter, int _verb, double *_weight) :
           myVerify(_Comm),
           callBLAS(),
           callFortran(),
           modalTool(_Comm),
           mySort(),
           MyComm(_Comm),
           K(KK),
           M(MM),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           normWeight(_weight),
           verbose(_verb),
           historyCount(0),
           resHistory(0),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           numResidual(0),
           orthoOp(0),
           outerIter(0),
           stifOp(0),
           timeMassOp(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timeResidual(0.0),
           timeStifOp(0.0)
           {

}


ModifiedARPACKm3::~ModifiedARPACKm3() {

  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

}


int ModifiedARPACKm3::solve(int numEigen, Epetra_MultiVector &Q, double *lambda) {

  // Computes the smallest eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = M X Lambda
  // 
  // using ModifiedARPACK (mode 3).
  //
  // The convergence test is performed outisde of ARPACK
  //
  //                      || Kx - Mx lambda || < tol*lambda
  //
  // The norm ||.|| is specified by the user through the array normWeight.
  //
  // Note that if M is not specified, then  K X = X Lambda is solved.
  // (using the mode for generalized eigenvalue problem).
  // 
  // Input variables:
  // 
  // numEigen  (integer) = Number of eigenmodes requested
  // 
  // Q (Epetra_MultiVector) = Initial search space
  //                   The number of columns of Q defines the size of search space (=NCV).
  //                   The rows of X are distributed across processors.
  //                   As a rule of thumb in ARPACK User's guide, NCV >= 2*numEigen.
  //                   At exit, the first numEigen locations contain the eigenvectors requested.
  // 
  // lambda (array of doubles) = Converged eigenvalues
  //                   The length of this array is equal to the number of columns in Q.
  //                   At exit, the first numEigen locations contain the eigenvalues requested.
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
  // info = -  8 >> numEigen must be smaller than the dimension of the matrix.
  //
  // info = - 30 >> MEMORY
  //
  // See ARPACK documentation for the meaning of INFO

  if (numEigen <= 0) {
    return 0;
  }

  int info = myVerify.inputArguments(numEigen, K, M, 0, Q, minimumSpaceDimension(numEigen));
  if (info < 0)
    return info;

  int myPid = MyComm.MyPID();

  int localSize = Q.MyLength();
  int NCV = Q.NumVectors();
  int knownEV = 0;

  if (NCV > Q.GlobalLength()) {
    if (numEigen >= Q.GlobalLength()) {
      cerr << endl;
      cerr << " !! The number of requested eigenvalues must be smaller than the dimension";
      cerr << " of the matrix !!\n";
      cerr << endl;
      return -8;
    }
    NCV = Q.GlobalLength();
  }

  // Get the weight for approximating the M-inverse norm
  Epetra_Vector *vectWeight = 0;
  if (normWeight) {
    vectWeight = new Epetra_Vector(View, Q.Map(), normWeight);
  }

  int localVerbose = verbose*(myPid == 0);

  // Define data for ARPACK
  //
  // UH (10/17/03) Note that workl is also used 
  //               * to store the eigenvectors of the tridiagonal matrix
  //               * as a workspace for DSTEQR
  //               * as a workspace for recovering the global eigenvectors

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  int ido = 0;

  int lwI = 22;
  int *wI = new (nothrow) int[lwI];
  if (wI == 0) {
    if (vectWeight)
      delete vectWeight;
    return -30;
  }
  memRequested += sizeof(int)*lwI/(1024.0*1024.0);

  int *iparam = wI;
  int *ipntr = wI + 11;

  int lworkl = NCV*(NCV+8);
  int lwD = lworkl + 4*localSize;
  double *wD = new (nothrow) double[lwD];
  if (wD == 0) {
    if (vectWeight)
      delete vectWeight;
    delete[] wI;
    return -30;
  }
  memRequested += sizeof(double)*(4*localSize+lworkl)/(1024.0*1024.0);

  double *pointer = wD;

  double *workl = pointer;
  pointer = pointer + lworkl;

  double *resid = pointer;
  pointer = pointer + localSize;

  double *workd = pointer;

  double *v = Q.Values();

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  iparam[1-1] = 1;
  iparam[3-1] = maxIterEigenSolve;
  iparam[7-1] = 3;

  // The fourth parameter forces to use the convergence test provided by ARPACK.
  // This requires a customization of ARPACK (provided by R. Lehoucq).

  iparam[4-1] = 1;

  Epetra_Vector v1(View, Q.Map(), workd);
  Epetra_Vector v2(View, Q.Map(), workd + localSize);
  Epetra_Vector v3(View, Q.Map(), workd + 2*localSize);

  // Define further storage for the new residual check
  // Use a block of vectors to compute the residuals more quickly.
  // Note that workd could be used if memory becomes an issue.
  int loopZ = (NCV > 10) ? 10 : NCV;

  int lwD2 = localSize + 2*NCV-1 + NCV;
  lwD2 += (M) ? 3*loopZ*localSize : 2*loopZ*localSize;
  double *wD2 = new (nothrow) double[lwD2];
  if (wD2 == 0) {
    if (vectWeight)
      delete vectWeight;
    delete[] wI;
    delete[] wD;
    return -30;
  }
  memRequested += sizeof(double)*lwD2/(1024.0*1024.0);

  pointer = wD2;
  
  // vTmp is used when ido = -1
  double *vTmp = pointer;
  pointer = pointer + localSize;

  // dd and ee stores the tridiagonal matrix.
  // Note that DSTEQR destroys the contents of the input arrays.
  double *dd = pointer;
  pointer = pointer + NCV;

  double *ee = pointer;
  pointer = pointer + NCV-1;

  double *vz = pointer;
  pointer = pointer + loopZ*localSize;
  Epetra_MultiVector approxEV(View, Q.Map(), vz, localSize, loopZ);

  double *kvz = pointer;
  pointer = pointer + loopZ*localSize;
  Epetra_MultiVector KapproxEV(View, Q.Map(), kvz, localSize, loopZ);

  double *mvz = (M) ? pointer : vz;
  pointer = (M) ? pointer + loopZ*localSize : pointer;
  Epetra_MultiVector MapproxEV(View, Q.Map(), mvz, localSize, loopZ);

  double *normR = pointer;

  // zz contains the eigenvectors of the tridiagonal matrix.
  // workt is a workspace for DSTEQR.
  // Note that zz and workt will use parts of workl.
  double *zz, *workt;

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  // Define an array to store the residuals history
  if (localVerbose > 1) {
    resHistory = new (nothrow) double[maxIterEigenSolve*NCV];
    if (resHistory == 0) {
      if (vectWeight)
        delete vectWeight;
      delete[] wI;
      delete[] wD;
      delete[] wD2;
      return -30;
    }
    historyCount = 0;
  }

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  if (localVerbose > 0) {
    cout << endl;
    cout << " *|* Problem: ";
    if (M) 
      cout << "K*Q = M*Q D ";
    else
      cout << "K*Q = Q D ";
    cout << endl;
    cout << " *|* Algorithm = ARPACK (Mode 3, modified such that user checks convergence)" << endl;
    cout << " *|* Number of requested eigenvalues = " << numEigen << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " *|* Tolerance for convergence = " << tolEigenSolve << endl;
    cout << " *|* Norm used for convergence: ";
    if (normWeight)
      cout << "weighted L2-norm with user-provided weights" << endl;
    else
      cout << "L^2-norm" << endl;
    cout << "\n -- Start iterations -- \n";
  }

#ifdef EPETRA_MPI
  Epetra_MpiComm *MPIComm = dynamic_cast<Epetra_MpiComm *>(const_cast<Epetra_Comm*>(&MyComm));
#endif

  timeOuterLoop -= MyWatch.WallTime();
  while (ido != 99) {

    highMem = (highMem > currentSize()) ? highMem : currentSize();

#ifdef EPETRA_MPI
    if (MPIComm)
      callFortran.PSAUPD(MPIComm->Comm(), &ido, 'G', localSize, "LM", numEigen, tolEigenSolve,
             resid, NCV, v, localSize, iparam, ipntr, workd, workl, lworkl, &info, 0);
    else
      callFortran.SAUPD(&ido, 'G', localSize, "LM", numEigen, tolEigenSolve, resid, NCV, v, 
             localSize, iparam, ipntr, workd, workl, lworkl, &info, 0);
#else
    callFortran.SAUPD(&ido, 'G', localSize, "LM", numEigen, tolEigenSolve, resid, NCV, v,
             localSize, iparam, ipntr, workd, workl, lworkl, &info, 0);
#endif

    if (ido == -1) {
      // Apply the mass matrix      
      v3.ResetView(workd + ipntr[0] - 1);
      v1.ResetView(vTmp);
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(v3, v1);
      else
        memcpy(v1.Values(), v3.Values(), localSize*sizeof(double));
      timeMassOp += MyWatch.WallTime();
      massOp += 1;
      // Solve the stiffness problem
      v2.ResetView(workd + ipntr[1] - 1);
      timeStifOp -= MyWatch.WallTime();
      K->ApplyInverse(v1, v2);
      timeStifOp += MyWatch.WallTime();
      stifOp += 1;
      continue;
    } // if (ido == -1)

    if (ido == 1) {
      // Solve the stiffness problem
      v1.ResetView(workd + ipntr[2] - 1);
      v2.ResetView(workd + ipntr[1] - 1);
      timeStifOp -= MyWatch.WallTime();
      K->ApplyInverse(v1, v2);
      timeStifOp += MyWatch.WallTime();
      stifOp += 1;
      continue;
    } // if (ido == 1)

    if (ido == 2) {
      // Apply the mass matrix      
      v1.ResetView(workd + ipntr[0] - 1);
      v2.ResetView(workd + ipntr[1] - 1);
      timeMassOp -= MyWatch.WallTime();
      if (M)
        M->Apply(v1, v2);
      else
        memcpy(v2.Values(), v1.Values(), localSize*sizeof(double));
      timeMassOp += MyWatch.WallTime();
      massOp += 1;
      continue;
    } // if (ido == 2)

    if (ido == 4) {
      timeResidual -= MyWatch.WallTime();
      // Copy the main diagonal of T
      memcpy(dd, workl + NCV + ipntr[4] - 1, NCV*sizeof(double));
      // Copy the lower diagonal of T
      memcpy(ee, workl + ipntr[4], (NCV-1)*sizeof(double));
      // Compute the eigenpairs of the tridiagonal matrix
      zz = workl + 4*NCV;
      workt = workl + 4*NCV + NCV*NCV;
      callFortran.STEQR('I', NCV, dd, ee, zz, NCV, workt, &info);
      if (info != 0) {
        if (localVerbose > 0) {
          cerr << endl;
          cerr << " Error with DSTEQR, info = " << info << endl;
          cerr << endl;
        }
        break;
      }
      // dd contains the eigenvalues in ascending order 
      // Check the residual of the proposed eigenvectors of (K, M)
      int ii, jz;
      iparam[4] = 0;
      for (jz = 0; jz < NCV; jz += loopZ) {
        int colZ = (jz + loopZ < NCV) ? loopZ : NCV - jz;
        callBLAS.GEMM('N', 'N', localSize, colZ, NCV, 1.0, v, localSize,
                      zz + jz*NCV, NCV, 0.0, vz, localSize);
        // Form the residuals
        if (M)
          M->Apply(approxEV, MapproxEV); 
        K->Apply(approxEV, KapproxEV); 
        for (ii = 0; ii < colZ; ++ii) {
          callBLAS.AXPY(localSize, -1.0/dd[ii+jz], MapproxEV.Values() + ii*localSize, 
                        KapproxEV.Values() + ii*localSize);
        }
        // Compute the norms of the residuals
        if (vectWeight) {
          KapproxEV.NormWeighted(*vectWeight, normR + jz);
        }
        else {
          KapproxEV.Norm2(normR + jz);
        }
        // Scale the norms of residuals with the eigenvalues
        for (ii = 0; ii < colZ; ++ii) {
          normR[ii+jz] = normR[ii+jz]*dd[ii+jz];
        }
        // Put the number of converged pairs in iparam[5-1]
        for (ii=0; ii<colZ; ++ii) {
          if (normR[ii+jz] < tolEigenSolve)
            iparam[4] += 1;
        }
      }
      timeResidual += MyWatch.WallTime();
      numResidual += NCV;
      outerIter += 1;
      if (localVerbose > 0) {
        cout << " Iteration " << outerIter;
        cout << " - Number of converged eigenvalues " << iparam[4] << endl;
      }
      if (localVerbose > 1) {
        memcpy(resHistory + historyCount, normR, NCV*sizeof(double));
        historyCount += NCV;
      }
      if (localVerbose > 3) {
        cout.precision(2);
        cout.setf(ios::scientific, ios::floatfield);
        for (ii=0; ii < NCV; ++ii) {
          cout << " Iteration " << outerIter;
          cout << " - Scaled Norm of Residual " << ii << " = " << normR[ii] << endl;
        }
        cout << endl;
        cout.precision(2);
        for (ii = 0; ii < NCV; ++ii) {
          cout << " Iteration " << outerIter << " - Ritz eigenvalue " << ii;
          cout.setf((fabs(dd[ii]) > 100) ? ios::scientific : ios::fixed, ios::floatfield);
          cout << " = " << 1.0/dd[ii] << endl;
        }
        cout << endl;
      }
    } // if (ido == 4)

  } // while (ido != 99)
  timeOuterLoop += MyWatch.WallTime();
  highMem = (highMem > currentSize()) ? highMem : currentSize();

  if (info < 0) {
    if (myPid == 0) {
      cerr << endl;
      cerr << " Error with DSAUPD, info = " << info << endl;
      cerr << endl;
    }
  }
  else {
    // Get the eigenvalues
    timePostProce -= MyWatch.WallTime();
    int ii, jj;
    double *pointer = workl + 4*NCV + NCV*NCV;
    for (ii=0; ii < localSize; ii += 3) {
      int nRow = (ii + 3 < localSize) ? 3 : localSize - ii;
      for (jj=0; jj<NCV; ++jj)
        memcpy(pointer + jj*nRow, v + ii + jj*localSize, nRow*sizeof(double));
      callBLAS.GEMM('N', 'N', nRow, NCV, NCV, 1.0, pointer, nRow, zz, NCV,
                    0.0, Q.Values() + ii, localSize);
    }
    // Put the converged eigenpairs at the beginning
    knownEV = 0;
    for (ii=0; ii < NCV; ++ii) {
      if (normR[ii] < tolEigenSolve) {
        lambda[knownEV] = 1.0/dd[ii];
        memcpy(Q.Values()+knownEV*localSize, Q.Values()+ii*localSize, localSize*sizeof(double));
        knownEV += 1;
        if (knownEV == Q.NumVectors())
          break;
      }
    }   
    // Sort the eigenpairs
    if (knownEV > 0) {
      mySort.sortScalars_Vectors(knownEV, lambda, Q.Values(), localSize);
    }
    timePostProce += MyWatch.WallTime();
  } // if (info < 0)

  if (info == 0) {
    orthoOp = iparam[11-1];
  }

  delete[] wI;
  delete[] wD;
  delete[] wD2;
  if (vectWeight)
    delete vectWeight;

  return (info == 0) ? knownEV : info;

}


void ModifiedARPACKm3::initializeCounters() {

  historyCount = 0;
  if (resHistory) {
    delete[] resHistory;
    resHistory = 0;
  }

  memRequested = 0.0;
  highMem = 0.0;

  massOp = 0;
  numResidual = 0;
  orthoOp = 0;
  outerIter = 0;
  stifOp = 0;

  timeMassOp = 0.0;
  timeOuterLoop = 0.0;
  timePostProce = 0.0;
  timeResidual = 0.0;
  timeStifOp = 0.0;

}


void ModifiedARPACKm3::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: ARPACK (Mode 3, modified such that user checks convergence)\n";
  }

}


void ModifiedARPACKm3::memoryInfo() const {

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
    cout << maxMemRequested << endl;
    cout << endl;
    cout << " High water mark in eigensolver                      = (EST) ";
    cout.width(6);
    cout << maxHighMem << endl;
    cout << endl;
  }

}


void ModifiedARPACKm3::historyInfo() const {

  if (resHistory) {
    int j;
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


void ModifiedARPACKm3::operationInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Operations ---\n\n";
    cout << " Total number of mass matrix multiplications      = ";
    cout.width(9);
    cout << massOp << endl;
    cout << " Total number of orthogonalizations               = ";
    cout.width(9);
    cout << orthoOp << endl;
    cout << " Total number of stiffness matrix operations      = ";
    cout.width(9);
    cout << stifOp << endl;
    cout << " Total number of computed eigen-residuals         = ";
    cout.width(9);
    cout << numResidual << endl;
    cout << endl;
    cout << " Total number of outer iterations                 = ";
    cout.width(9);
    cout << outerIter << endl;
    cout << endl;
  }

}


void ModifiedARPACKm3::timeInfo() const {

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
    cout << "       Time for stiffness matrix solves          = ";
    cout.width(9);
    cout << timeStifOp << " s     ";
    cout.width(5);
    cout << 100*timeStifOp/timeOuterLoop << " %\n";
    cout << "       Time for residual computations            = ";
    cout.width(9);
    cout << timeResidual << " s     ";
    cout.width(5);
    cout << 100*timeResidual/timeOuterLoop << " %\n";
    cout << endl;
    cout << " Total time for post-processing                  = ";
    cout.width(9);
    cout << timePostProce << " s\n";
    cout << endl;
  } // if (myId == 0)

}


