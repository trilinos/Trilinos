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

#include "ARPACKm1.h"


ARPACKm1::ARPACKm1(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                     double _tol, int _maxIter, int _verb) :
           myVerify(_Comm),
           callFortran(),
           MyComm(_Comm),
           K(KK),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           which("LM"),
           verbose(_verb),
           memRequested(0.0),
           highMem(0.0),
           orthoOp(0),
           outerIter(0),
           stifOp(0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timeStifOp(0.0)
           {

}


ARPACKm1::ARPACKm1(const Epetra_Comm &_Comm, const Epetra_Operator *KK, char *_which,
                   double _tol, int _maxIter, int _verb) :
           myVerify(_Comm),
           callFortran(),
           MyComm(_Comm),
           K(KK),
           MyWatch(_Comm),
           tolEigenSolve(_tol),
           maxIterEigenSolve(_maxIter),
           which(_which),
           verbose(_verb),
           memRequested(0.0),
           highMem(0.0),
           orthoOp(0),
           outerIter(0),
           stifOp(0),
           timeOuterLoop(0.0),
           timePostProce(0.0),
           timeStifOp(0.0)
           {

}


int ARPACKm1::solve(int numEigen, Epetra_MultiVector &Q, double *lambda) {

  return ARPACKm1::reSolve(numEigen, Q, lambda, 0);

}


int ARPACKm1::reSolve(int numEigen, Epetra_MultiVector &Q, double *lambda, int startingEV) {

  // Computes eigenvalues and the corresponding eigenvectors
  // of the generalized eigenvalue problem
  // 
  //      K X = X Lambda
  // 
  // using ARPACK (mode 1).
  //
  // The convergence test is provided by ARPACK.
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
  // startingEV (integer) = Number of eigenmodes already stored in Q
  //                   A linear combination of these vectors is made to define the starting
  //                   vector, placed in resid.
  //
  // Return information on status of computation
  // 
  // info >=   0 >> Number of converged eigenpairs at the end of computation
  // 
  // // Failure due to input arguments
  // 
  // info = -  1 >> The stiffness matrix K has not been specified.
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

  if (numEigen <= startingEV) {
    return numEigen;
  }

  int info = myVerify.inputArguments(numEigen, K, 0, 0, Q, minimumSpaceDimension(numEigen));
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

  if (startingEV > 0) {
    // Define the initial starting vector
    memset(resid, 0, localSize*sizeof(double));
    for (int jj = 0; jj < startingEV; ++jj)
      for (int ii = 0; ii < localSize; ++ii)
         resid[ii] += v[ii + jj*localSize];
    info = 1;
  }

  iparam[1-1] = 1;
  iparam[3-1] = maxIterEigenSolve;
  iparam[7-1] = 1;

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
    cout << "K*Q = Q D ";
    cout << endl;
    cout << " *|* Algorithm = ARPACK (mode 1)" << endl;
    cout << " *|* Number of requested eigenvalues = " << numEigen << endl;
    cout.precision(2);
    cout.setf(ios::scientific, ios::floatfield);
    cout << " *|* Tolerance for convergence = " << tolEigenSolve << endl;
    if (startingEV > 0)
      cout << " *|* User-defined starting vector (Combination of " << startingEV << " vectors)\n";
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
      callFortran.PSAUPD(MPIComm->Comm(), &ido, 'I', localSize, which, numEigen, tolEigenSolve,
             resid, NCV, v, localSize, iparam, ipntr, workd, workl, lworkl, &info, localVerbose);
    else
      callFortran.SAUPD(&ido, 'I', localSize, which, numEigen, tolEigenSolve, resid, NCV, v, 
             localSize, iparam, ipntr, workd, workl, lworkl, &info, localVerbose);
#else
    callFortran.SAUPD(&ido, 'I', localSize,  which, numEigen, tolEigenSolve, resid, NCV, v,
             localSize, iparam, ipntr, workd, workl, lworkl, &info, localVerbose);
#endif

    if ((ido == 1) || (ido == -1)) {
      // Apply the matrix
      v1.ResetView(workd + ipntr[0] - 1);
      v2.ResetView(workd + ipntr[1] - 1);
      timeStifOp -= MyWatch.WallTime();
      K->Apply(v1, v2);
      timeStifOp += MyWatch.WallTime();
      stifOp += 1;
      continue;
    } // if ((ido == 1) || (ido == -1))

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
      callFortran.PSEUPD(MPIComm->Comm(), 1, 'A', select, lambda, v, localSize, sigma, 'I',
            localSize, which, numEigen, tolEigenSolve, resid, NCV, v, localSize, iparam, ipntr, 
            workd, workl, lworkl, &info);
    else
      callFortran.SEUPD(1, 'A', select, lambda, v, localSize, sigma, 'I', localSize, which,
            numEigen, tolEigenSolve, resid, NCV, v, localSize, iparam, ipntr, workd, workl,
            lworkl, &info);
#else
    callFortran.SEUPD(1, 'A', select, lambda, v, localSize, sigma, 'I', localSize, which,
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


void ARPACKm1::initializeCounters() {

  memRequested = 0.0;
  highMem = 0.0;

  orthoOp = 0;
  outerIter = 0;
  stifOp = 0;

  timeOuterLoop = 0.0;
  timePostProce = 0.0;
  timeStifOp = 0.0;

}


void ARPACKm1::algorithmInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid ==0) {
    cout << " Algorithm: ARPACK (Mode 1)\n";
  }

}


void ARPACKm1::memoryInfo() const {

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


void ARPACKm1::operationInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Operations ---\n\n";
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


void ARPACKm1::timeInfo() const {

  int myPid = MyComm.MyPID();
  
  if (myPid == 0) {
    cout << " --- Timings ---\n\n";
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(2);
    cout << " Total time for outer loop                       = ";
    cout.width(9);
    cout << timeOuterLoop << " s" << endl;
    cout << "       Time for stiffness matrix                 = ";
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


