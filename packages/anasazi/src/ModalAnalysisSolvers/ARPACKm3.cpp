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

#include "ARPACKm3.h"


ARPACKm3::ARPACKm3(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     double _tol, int _maxIter, int _verb) :
           myVerify(_Comm),
           callFortran(),
           MyComm(_Comm),
           K(KK),
           M(0),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           verbose(_verb),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           orthoOp(0),
           outerIter(0),
           stifOp(0),
           timeMassOp(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timeStifOp(0.0)
           {

}


ARPACKm3::ARPACKm3(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     const Epetra_Operator *MM,
                     double _tol, int _maxIter, int _verb) :
           myVerify(_Comm),
           callFortran(),
           MyComm(_Comm),
           K(KK),
           M(MM),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           verbose(_verb),
           memRequested(0.0),
           highMem(0.0),
           massOp(0),
           orthoOp(0),
           outerIter(0),
           stifOp(0),
           timeMassOp(0.0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timeStifOp(0.0)
           {

}


int ARPACKm3::solve(int numEigen, Epetra_MultiVector &Q, double *lambda) {

  // Computes the smallest eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = M X Lambda
  // 
  // using ARPACK (mode 3).
  //
  // The convergence test is provided by ARPACK.
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

  int localVerbose = verbose*(myPid == 0);

  // Define data for ARPACK
  highMem = (highMem > currentSize()) ? highMem : currentSize();

  int ido = 0;

  int lwI = 22 + NCV;
  int *wI = new (nothrow) int[lwI];
  if (wI == 0) {
    return -30;
  }
  memRequested += sizeof(int)*lwI/(1024.0*1024.0);

  int *iparam = wI;
  int *ipntr = wI + 11;
  int *select = wI + 22;

  int lworkl = NCV*(NCV+8);
  int lwD = lworkl + 4*localSize;
  double *wD = new (nothrow) double[lwD];
  if (wD == 0) {
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

  double sigma = 0.0;

  iparam[1-1] = 1;
  iparam[3-1] = maxIterEigenSolve;
  iparam[7-1] = 3;

  // The fourth parameter forces to use the convergence test provided by ARPACK.
  // This requires a customization of ARPACK (provided by R. Lehoucq).

  iparam[4-1] = 0;

  Epetra_Vector v1(View, Q.Map(), workd);
  Epetra_Vector v2(View, Q.Map(), workd + localSize);
  Epetra_Vector v3(View, Q.Map(), workd + 2*localSize);

  double *vTmp = new (nothrow) double[localSize];
  if (vTmp == 0) {
    delete[] wI;
    delete[] wD;
    return -30;
  }
  memRequested += sizeof(double)*localSize/(1024.0*1024.0);

  highMem = (highMem > currentSize()) ? highMem : currentSize();

  if (localVerbose > 0) {
    cout << endl;
    cout << " *|* Problem: ";
    if (M) 
      cout << "K*Q = M*Q D ";
    else
      cout << "K*Q = Q D ";
    cout << endl;
    cout << " *|* Algorithm = ARPACK (mode 3)" << endl;
    cout << " *|* Number of requested eigenvalues = " << numEigen << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " *|* Tolerance for convergence = " << tolEigenSolve << endl;
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
             resid, NCV, v, localSize, iparam, ipntr, workd, workl, lworkl, &info, localVerbose);
    else
      callFortran.SAUPD(&ido, 'G', localSize, "LM", numEigen, tolEigenSolve, resid, NCV, v, 
             localSize, iparam, ipntr, workd, workl, lworkl, &info, localVerbose);
#else
    callFortran.SAUPD(&ido, 'G', localSize, "LM", numEigen, tolEigenSolve, resid, NCV, v,
             localSize, iparam, ipntr, workd, workl, lworkl, &info, localVerbose);
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

    // Compute the eigenvectors
    timePostProce -= MyWatch.WallTime();
#ifdef EPETRA_MPI
    if (MPIComm)
      callFortran.PSEUPD(MPIComm->Comm(), 1, 'A', select, lambda, v, localSize, sigma, 'G',
            localSize, "LM", numEigen, tolEigenSolve, resid, NCV, v, localSize, iparam, ipntr, 
            workd, workl, lworkl, &info);
    else
      callFortran.SEUPD(1, 'A', select, lambda, v, localSize, sigma, 'G', localSize, "LM",
            numEigen, tolEigenSolve, resid, NCV, v, localSize, iparam, ipntr, workd, workl,
            lworkl, &info);
#else
    callFortran.SEUPD(1, 'A', select, lambda, v, localSize, sigma, 'G', localSize, "LM",
          numEigen, tolEigenSolve, resid, NCV, v, localSize, iparam, ipntr, workd, workl,
          lworkl, &info);
#endif
    timePostProce += MyWatch.WallTime();
    highMem = (highMem > currentSize()) ? highMem : currentSize();

    // Treat the error
    if (info != 0) {
      if (myPid == 0) {
        cerr << endl;
        cerr << " Error with DSEUPD, info = " << info << endl;
        cerr << endl;
      }
    }

  } // if (info < 0)

  if (info == 0) {
    outerIter = iparam[3-1];
    knownEV = iparam[5-1];
    orthoOp = iparam[11-1];
  }

  delete[] wI;
  delete[] wD;
  delete[] vTmp;

  return (info == 0) ? knownEV : info;

}


void ARPACKm3::initializeCounters() {

  memRequested = 0.0;
  highMem = 0.0;

  massOp = 0;
  orthoOp = 0;
  outerIter = 0;
  stifOp = 0;

  timeMassOp = 0.0;
  timeOuterLoop = 0.0;
  timePostProce = 0.0;
  timeStifOp = 0.0;

}


void ARPACKm3::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: ARPACK (Mode 3)\n";
  }

}


void ARPACKm3::memoryInfo() const {

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


void ARPACKm3::operationInfo() const {

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
    cout << endl;
    cout << " Total number of outer iterations                 = ";
    cout.width(9);
    cout << outerIter << endl;
    cout << endl;
  }

}


void ARPACKm3::timeInfo() const {

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
    cout << endl;
    cout << " Total time for post-processing                  = ";
    cout.width(9);
    cout << timePostProce << " s\n";
    cout << endl;
  } // if (myId == 0)

}


