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

#include "ModalTools.h"


ModalTools::ModalTools(const Epetra_Comm &_Comm) :
         callFortran(),
         callBLAS(),
         callLAPACK(),
         MyComm(_Comm),
         MyWatch(_Comm),
         eps(0.0),
         timeQtMult(0.0),
         timeQMult(0.0),
         timeProj_MassMult(0.0),
         timeNorm_MassMult(0.0),
         timeProj(0.0),
         timeNorm(0.0),
         numProj_MassMult(0),
         numNorm_MassMult(0) {

  callLAPACK.LAMCH('E', eps);

}


int ModalTools::makeSimpleLumpedMass(const Epetra_Operator *M, double *normWeight) const {

  // Return, in the array normWeight, weights value such that
  // the function Epetra_MultiVector::NormWeighted computes
  //
  //              || . || = ( N / ( 1^T M 1) )^{1/2} * || . ||_{2}
  //
  // where 1 is the vector with unit entries and N is the size of the problem
  //
  // normWeight : Array of double (length = # of unknowns on this processor)
  //              Contains at exit the weights
  //
  // Note: When M is 0, the function does not do anything

  if (M == 0)
    return -1;

  int row = (M->OperatorDomainMap()).NumMyElements();

  double *zVal = new double[row];

  Epetra_Vector z(View, M->OperatorDomainMap(), zVal);
  Epetra_Vector Mz(View, M->OperatorDomainMap(), normWeight);

  z.PutScalar(1.0);
  M->Apply(z, Mz);

  delete[] zVal;

  int i;
  double rho = 0.0;
  for (i = 0; i < row; ++i)
    rho += normWeight[i];

  normWeight[0] = rho;
  MyComm.SumAll(normWeight, &rho, 1);

  int info = 0;
  if (rho <= 0.0) {
    info = -2;
    if (MyComm.MyPID() == 0)
      cerr << "\n !!! The mass matrix gives a negative mass: " << rho << " !!! \n\n";
    return info;
  }
  rho = rho/(M->OperatorDomainMap()).NumGlobalElements();

  // 
  // Note that the definition of the weighted norm in Epetra forces this weight definition
  // UH 09/03/03
  //

  rho = sqrt(rho/(M->OperatorDomainMap()).NumGlobalElements());

  for (i = 0; i < row; ++i)
   normWeight[i] = rho;

  return info;

}


int ModalTools::massOrthonormalize(Epetra_MultiVector &X, Epetra_MultiVector &MX, 
                                const Epetra_Operator *M, const Epetra_MultiVector &Q,
                                int howMany, int type, double *WS, double kappa) {

  // For the inner product defined by the operator M or the identity (M = 0)
  //   -> Orthogonalize X against Q
  //   -> Orthonormalize X 
  // Modify MX accordingly
  // WS is used as a workspace (size: (# of columns in X)*(# of rows in X))
  //
  // Note that when M is 0, MX is not referenced
  //
  // Parameter variables
  //
  // X  : Vectors to be transformed
  //
  // MX : Image of the block vector X by the mass matrix
  //
  // Q  : Vectors to orthogonalize against
  //
  // howMany : Number of vectors of X to orthogonalized
  //           If this number is smaller than the total number of vectors in X,
  //           then it is assumed that the last "howMany" vectors are not orthogonal
  //           while the other vectors in X are othogonal to Q and orthonormal.
  //
  // type = 0 (default) > Performs both operations
  // type = 1           > Performs Q^T M X = 0
  // type = 2           > Performs X^T M X = I
  //
  // WS   = Working space (default value = 0)
  //
  // kappa= Coefficient determining when to perform a second Gram-Schmidt step
  //        Default value = 1.5625 = (1.25)^2 (as suggested in Parlett's book)
  //
  // Output the status of the computation
  //
  // info =   0 >> Success.
  //
  // info >   0 >> Indicate how many vectors have been tried to avoid rank deficiency for X 
  //
  // info =  -1 >> Failure >> X has zero columns
  //                       >> It happens when # col of X     > # rows of X 
  //                       >> It happens when # col of [Q X] > # rows of X 
  //                       >> It happens when no good random vectors could be found

  int info = 0;

  // Orthogonalize X against Q
  timeProj -= MyWatch.WallTime();
  if (type != 2) {

    int xc = howMany;
    int xr = X.MyLength();
    int qc = Q.NumVectors();

    Epetra_MultiVector XX(View, X, X.NumVectors()-howMany, howMany);

    Epetra_MultiVector MXX(View, X.Map(), (M) ? MX.Values() + xr*(MX.NumVectors()-howMany)
                                              : X.Values() + xr*(X.NumVectors()-howMany),
                           xr, howMany);

    // Perform the Gram-Schmidt transformation for a block of vectors

    // Compute the initial M-norms
    double *oldDot = new double[xc];
    MXX.Dot(XX, oldDot);

    // Define the product Q^T * (M*X)
    double *qTmx = new double[2*qc*xc];

    // Multiply Q' with MX
    timeQtMult -= MyWatch.WallTime();
    callBLAS.GEMM('T', 'N', qc, xc, xr, 1.0, Q.Values(), xr, MXX.Values(), xr, 
                 0.0, qTmx + qc*xc, qc);
    MyComm.SumAll(qTmx + qc*xc, qTmx, qc*xc);
    timeQtMult += MyWatch.WallTime();

    // Multiply by Q and substract the result in X
    timeQMult -= MyWatch.WallTime();
    callBLAS.GEMM('N', 'N', xr, xc, qc, -1.0, Q.Values(), xr, qTmx, qc, 
                  1.0, XX.Values(), xr); 
    timeQMult += MyWatch.WallTime();

    // Update MX
    if (M) {
      if ((qc >= xc) || (WS == 0)) {
        timeProj_MassMult -= MyWatch.WallTime();
        M->Apply(XX, MXX);
        timeProj_MassMult += MyWatch.WallTime();
        numProj_MassMult += xc;
      }
      else {
        Epetra_MultiVector MQ(View, Q.Map(), WS, Q.MyLength(), qc);
        timeProj_MassMult -= MyWatch.WallTime();
        M->Apply(Q, MQ);
        timeProj_MassMult += MyWatch.WallTime();
        numProj_MassMult += qc;
        callBLAS.GEMM('N', 'N', xr, xc, qc, -1.0, MQ.Values(), xr, qTmx, qc, 
                      1.0, MXX.Values(), xr); 
      }  // if ((qc >= xc) || (WS == 0))
    } // if (M)

    double newDot = 0.0;
    int j;
    for (j = 0; j < xc; ++j) {

      MXX(j)->Dot(*(XX(j)), &newDot);

      if (kappa*newDot < oldDot[j]) {

        // Apply another step of classical Gram-Schmidt
        timeQtMult -= MyWatch.WallTime();
        callBLAS.GEMM('T', 'N', qc, xc, xr, 1.0, Q.Values(), xr, MXX.Values(), xr, 
                      0.0, qTmx + qc*xc, qc);
        MyComm.SumAll(qTmx + qc*xc, qTmx, qc*xc);
        timeQtMult += MyWatch.WallTime();

        timeQMult -= MyWatch.WallTime();
        callBLAS.GEMM('N', 'N', xr, xc, qc, -1.0, Q.Values(), xr, qTmx, qc, 
                      1.0, XX.Values(), xr); 
        timeQMult += MyWatch.WallTime();

        // Update MX
        if (M) {
          if ((qc >= xc) || (WS == 0)) {
            timeProj_MassMult -= MyWatch.WallTime();
            M->Apply(XX, MXX);
            timeProj_MassMult += MyWatch.WallTime();
            numProj_MassMult += xc;
          }
          else {
            Epetra_MultiVector MQ(View, Q.Map(), WS, Q.MyLength(), qc);
            timeProj_MassMult -= MyWatch.WallTime();
            M->Apply(Q, MQ);
            timeProj_MassMult += MyWatch.WallTime();
            numProj_MassMult += qc;
            callBLAS.GEMM('N', 'N', xr, xc, qc, -1.0, MQ.Values(), xr, qTmx, qc, 
                          1.0, MXX.Values(), xr); 
          } // if ((qc >= xc) || (WS == 0))
        } // if (M)

        break;
      } // if (kappa*newDot < oldDot[j])
    } // for (j = 0; j < xc; ++j)

    delete[] qTmx;
    delete[] oldDot;

  } // if (type != 2)
  timeProj += MyWatch.WallTime();

  // Orthonormalize X 
  timeNorm -= MyWatch.WallTime();
  if (type != 1) {

    int j;
    int xc = X.NumVectors();
    int xr = X.MyLength();
    int globalSize = X.GlobalLength();
    int shift = (type == 2) ? 0 : Q.NumVectors();
    int mxc = (M) ? MX.NumVectors() : X.NumVectors();

    bool allocated = false;
    if (WS == 0) {
      allocated = true;
      WS = new double[xr];
    }

    double *oldMXj = WS;
    double *MXX = (M) ? MX.Values() : X.Values();
    double *product = new double[2*xc];

    double dTmp;

    for (j = 0; j < howMany; ++j) {

      int numX = xc - howMany + j;
      int numMX = mxc - howMany + j;

      // Put zero vectors in X when we are exceeding the space dimension
      if (numX + shift >= globalSize) {
        Epetra_Vector XXj(View, X, numX);
        XXj.PutScalar(0.0);
        if (M) {
          Epetra_Vector MXXj(View, MX, numMX);
          MXXj.PutScalar(0.0);
        }
        info = -1;
      }

      int numTrials;
      bool rankDef = true;
      for (numTrials = 0; numTrials < 10; ++numTrials) {

        double *Xj = X.Values() + xr*numX;
        double *MXj = MXX + xr*numMX;

        double oldDot = 0.0;
        dTmp = callBLAS.DOT(xr, Xj, MXj);
        MyComm.SumAll(&dTmp, &oldDot, 1);
      
        memcpy(oldMXj, MXj, xr*sizeof(double));

        if (numX > 0) {

          // Apply the first Gram-Schmidt

          callBLAS.GEMV('T', xr, numX, 1.0, X.Values(), xr, MXj, 0.0, product + xc);
          MyComm.SumAll(product + xc, product, numX);
          callBLAS.GEMV('N', xr, numX, -1.0, X.Values(), xr, product, 1.0, Xj);
          if (M) {
            if (xc == mxc) {
              callBLAS.GEMV('N', xr, numX, -1.0, MXX, xr, product, 1.0, MXj);
            }
            else {
              Epetra_Vector XXj(View, X, numX);
              Epetra_Vector MXXj(View, MX, numMX);
              timeNorm_MassMult -= MyWatch.WallTime();
              M->Apply(XXj, MXXj);
              timeNorm_MassMult += MyWatch.WallTime();
              numNorm_MassMult += 1;
            }
          }

          double dot = 0.0;
          dTmp = callBLAS.DOT(xr, Xj, MXj);
          MyComm.SumAll(&dTmp, &dot, 1);

          if (kappa*dot < oldDot) {
            callBLAS.GEMV('T', xr, numX, 1.0, X.Values(), xr, MXj, 0.0, product + xc);
            MyComm.SumAll(product + xc, product, numX);
            callBLAS.GEMV('N', xr, numX, -1.0, X.Values(), xr, product, 1.0, Xj);
            if (M) {
              if (xc == mxc) {
                callBLAS.GEMV('N', xr, numX, -1.0, MXX, xr, product, 1.0, MXj);
              }
              else {
                Epetra_Vector XXj(View, X, numX);
                Epetra_Vector MXXj(View, MX, numMX);
                timeNorm_MassMult -= MyWatch.WallTime();
                M->Apply(XXj, MXXj);
                timeNorm_MassMult += MyWatch.WallTime();
                numNorm_MassMult += 1;
              }
            }
          } // if (kappa*dot < oldDot)

        } // if (numX > 0)

        double norm = 0.0;
        dTmp = callBLAS.DOT(xr, Xj, oldMXj);
        MyComm.SumAll(&dTmp, &norm, 1);

        if (norm > oldDot*eps*eps) {
          norm = 1.0/sqrt(norm);
          callBLAS.SCAL(xr, norm, Xj);
          if (M)
           callBLAS.SCAL(xr, norm, MXj);
          rankDef = false;
          break;
        }
        else {
          info += 1;
          Epetra_Vector XXj(View, X, numX);
          XXj.Random();
          Epetra_Vector MXXj(View, MX, numMX);
          if (M) {
            timeNorm_MassMult -= MyWatch.WallTime();
            M->Apply(XXj, MXXj);
            timeNorm_MassMult += MyWatch.WallTime();
            numNorm_MassMult += 1;
          }
          if (type == 0)
            massOrthonormalize(XXj, MXXj, M, Q, 1, 1, WS, kappa);
        } // if (norm > oldDot*eps*eps)

      }  // for (numTrials = 0; numTrials < 10; ++numTrials)
  
      if (rankDef == true) {
        Epetra_Vector XXj(View, X, numX);
        XXj.PutScalar(0.0);
        if (M) {
          Epetra_Vector MXXj(View, MX, numMX);
          MXXj.PutScalar(0.0);
        }
        info = -1;
        break;
      }
  
    } // for (j = 0; j < howMany; ++j)

    delete[] product;

    if (allocated == true) {
      delete[] WS;
    }

  } // if (type != 1)
  timeNorm += MyWatch.WallTime();

  return info;

}


void ModalTools::localProjection(int numRow, int numCol, int dotLength,
                         double *U, int ldU, double *MatV, int ldV, 
                         double *UtMatV, int ldUtMatV, double *work) const {

  // This routine forms a projected matrix of a matrix Mat onto the subspace
  // spanned by U and V
  //
  // numRow = Number of columns in U (or number of rows in U^T) (input)
  //
  // numCol = Number of columns in V (input)
  //
  // dotLength = Local length of vectors U and V (input)
  //
  // U = Array of double storing the vectors U (input)
  //
  // ldU = Leading dimension in U (input)
  //
  // MatV = Array of double storing the vectors Mat*V (input)
  //
  // ldMatV = Leading dimension in MatV (input)
  //
  // UtMatV = Array of double storing the projected matrix U^T * Mat * V (output)
  //
  // ldUtMatV = Leading dimension in UtMatV (input)
  //
  // work = Workspace (size >= 2*numRow*numCol)

  int j;
  int offSet = numRow*numCol;

  callBLAS.GEMM('T', 'N', numRow, numCol, dotLength, 1.0, U, ldU, MatV, ldV, 
       0.0, work + offSet, numRow);
  MyComm.SumAll(work + offSet, work, offSet);
  
  double *source = work;
  double *result = UtMatV;
  int howMany =  numRow*sizeof(double);
  
  for (j = 0; j < numCol; ++j) {
    memcpy(result, source, howMany);
    // Shift the pointers
    source = source + numRow;
    result = result + ldUtMatV;
  }

}


int ModalTools::directSolver(int size, double *KK, int ldK, double *MM, int ldM,
             int &nev, double *EV, int ldEV, double *theta, int verbose, int esType) const {

  // Routine for computing the first NEV generalized eigenpairs of the symmetric pencil (KK, MM)
  //
  // Parameter variables:
  //
  // size : Size of the matrices KK and MM
  //
  // KK : Symmetric "stiffness" matrix (Size = size x ldK)
  //
  // ldK : Leading dimension of KK (ldK >= size)
  //
  // MM : Symmetric Positive "mass" matrix  (Size = size x ldM)
  //
  // ldM : Leading dimension of MM (ldM >= size)
  //
  // nev : Number of the smallest eigenvalues requested (input)
  //       Number of the smallest computed eigenvalues (output)
  //
  // EV : Array to store the eigenvectors (Size = nev x ldEV)
  //
  // ldEV : Leading dimension of EV (ldEV >= size)
  //
  // theta : Array to store the eigenvalues (Size = nev)
  // 
  // verbose : Flag to print information on the computation
  //
  // esType : Flag to select the algorithm
  //
  // esType =  0 (default) Uses LAPACK routine (Cholesky factorization of MM)
  //                       with deflation of MM to get orthonormality of 
  //                       eigenvectors (S^T MM S = I)
  //
  // esType =  1           Uses LAPACK routine (Cholesky factorization of MM)
  //                       (no check of orthonormality)
  //
  // esType = 10           Uses LAPACK routine for simple eigenproblem on KK
  //                       (MM is not referenced in this case)
  //
  // Note: The code accesses only the upper triangular part of KK and MM.
  //
  // Return the integer info on the status of the computation
  //
  // info = 0 >> Success
  //
  // info = - 20 >> Failure in LAPACK routine

  // Define local arrays

  double *kp = 0;
  double *mp = 0;
  double *tt = 0;

  double *U = 0;

  int i, j;
  int rank = 0;
  int info = 0;
  int tmpI;

  int NB = 5 + callFortran.LAENV(1, "dsytrd", "u", size, -1, -1, -1, 6, 1);
  int lwork = size*NB;
  double *work = new double[lwork];

//  double tol = sqrt(eps);
  double tol = 1e-12;

  switch (esType) {

      default:
      case 0:

      // Use the Cholesky factorization of MM to compute the generalized eigenvectors

      // Define storage
      kp = new double[size*size];
      mp = new double[size*size];
      tt = new double[size];
      U = new double[size*size];

      if (verbose > 4)
        cout << endl;

      // Copy KK & MM
      tmpI = sizeof(double);
      for (rank = size; rank > 0; --rank) {
        memset(kp, 0, size*size*tmpI);
        for (i = 0; i < rank; ++i) {
          memcpy(kp + i*size, KK + i*ldK, (i+1)*tmpI);
          memcpy(mp + i*size, MM + i*ldM, (i+1)*tmpI);
        }
        // Solve the generalized eigenproblem with LAPACK
        info = 0;
        callFortran.SYGV(1, 'V', 'U', rank, kp, size, mp, size, tt, work, lwork, &info);
        // Treat error messages
        if (info < 0) {
          if (verbose > 0) {
            cerr << endl;
            cerr << " In DSYGV, argument " << -info << "has an illegal value.\n";
            cerr << endl;
          }
          return -20;
        }
        if (info > 0) {
          if (info > rank)
            rank = info - rank;
          continue;
        }
        // Check the quality of eigenvectors
        for (i = 0; i < size; ++i) {
          memcpy(mp + i*size, MM + i*ldM, (i+1)*tmpI);
          for (j = 0; j < i; ++j)
            mp[i + j*size] = mp[j + i*size];
        }
        callBLAS.GEMM('N', 'N', size, rank, size, 1.0, mp, size, kp, size, 0.0, U, size);
        callBLAS.GEMM('T', 'N', rank, rank, size, 1.0, kp, size, U, size, 0.0, mp, rank);
        double maxNorm = 0.0;
        double maxOrth = 0.0;
        for (i = 0; i < rank; ++i) {
          for (j = i; j < rank; ++j) {
            if (j == i) {
              maxNorm = (fabs(mp[j+i*rank]-1.0) > maxNorm) ? fabs(mp[j+i*rank]-1.0) : maxNorm;
            }
            else {
              maxOrth = (fabs(mp[j+i*rank]) > maxOrth) ? fabs(mp[j+i*rank]) : maxOrth;
            }
          }
        }
        if (verbose > 4) {
          cout << " >> Local eigensolve >> Size: " << rank;
          cout.precision(2);
          cout.setf(ios::scientific, ios::floatfield);
          cout << " Normalization error: " << maxNorm;
          cout << " Orthogonality error: " << maxOrth;
          cout << endl;
        }
        if ((maxNorm <= tol) && (maxOrth <= tol))
          break;
      } // for (rank = size; rank > 0; --rank)

      if (verbose > 4)
        cout << endl;

      // Copy the eigenvectors
      memset(EV, 0, nev*ldEV*tmpI);
      nev = (rank < nev) ? rank : nev;
      memcpy(theta, tt, nev*tmpI);
      tmpI = rank*tmpI;
      for (i = 0; i < nev; ++i) {
        memcpy(EV + i*ldEV, kp + i*size, tmpI);
      }

      break;

      case 1:

      // Use the Cholesky factorization of MM to compute the generalized eigenvectors

      // Define storage
      kp = new double[size*size];
      mp = new double[size*size];
      tt = new double[size];

      // Copy KK & MM
      tmpI = sizeof(double);
      for (i = 0; i < size; ++i) {
        memcpy(kp + i*size, KK + i*ldK, (i+1)*tmpI);
        memcpy(mp + i*size, MM + i*ldM, (i+1)*tmpI);
      }

      // Solve the generalized eigenproblem with LAPACK
      info = 0;
      callFortran.SYGV(1, 'V', 'U', size, kp, size, mp, size, tt, work, lwork, &info);

      // Treat error messages
      if (info < 0) {
        if (verbose > 0) {
          cerr << endl;
          cerr << " In DSYGV, argument " << -info << "has an illegal value.\n";
          cerr << endl;
        }
        return -20;
      }
      if (info > 0) {
        if (info > size)
          nev = 0;
        else {
          if (verbose > 0) {
            cerr << endl;
            cerr << " In DSYGV, DPOTRF or DSYEV returned an error code (" << info << ").\n";
            cerr << endl;
          }
          return -20; 
        }
      }

      // Copy the eigenvectors
      memcpy(theta, tt, nev*tmpI);
      tmpI = size*tmpI;
      for (i = 0; i < nev; ++i) {
        memcpy(EV + i*ldEV, kp + i*size, tmpI);
      }

      break;

      case 10:

      // Simple eigenproblem

      // Define storage
      kp = new double[size*size];
      tt = new double[size];

      // Copy KK
      tmpI = sizeof(double);
      for (i=0; i<size; ++i) {
        memcpy(kp + i*size, KK + i*ldK, (i+1)*tmpI);
      }

      // Solve the generalized eigenproblem with LAPACK
      callFortran.SYEV('V', 'U', size, kp, size, tt, work, lwork, &info);

      // Treat error messages
      if (info != 0) {
        if (verbose > 0) {
          cerr << endl;
          if (info < 0) 
            cerr << " In DSYEV, argument " << -info << " has an illegal value\n";
          else
            cerr << " In DSYEV, the algorithm failed to converge (" << info << ").\n";
          cerr << endl;
        }
        info = -20;
        break;
      }

      // Copy the eigenvectors
      memcpy(theta, tt, nev*tmpI);
      tmpI = size*tmpI;
      for (i = 0; i < nev; ++i) {
        memcpy(EV + i*ldEV, kp + i*size, tmpI);
      }

      break;

  }

  // Clear the memory

  if (kp)
    delete[] kp;
  if (mp)
    delete[] mp;
  if (tt)
    delete[] tt;
  if (work)
    delete[] work;
  if (U)
    delete[] U;

  return info;

} 


