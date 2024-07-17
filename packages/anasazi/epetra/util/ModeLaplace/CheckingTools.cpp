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

#include "AnasaziConfigDefs.hpp"
#include "CheckingTools.h"


CheckingTools::CheckingTools(const Epetra_Comm &_Comm) :
         MyComm(_Comm) {

}


double CheckingTools::errorOrthogonality(const Epetra_MultiVector *X, 
                      const Epetra_MultiVector *R, const Epetra_Operator *M) const {

  // Return the maximum value of R_i^T * M * X_j / || MR_i || || X_j ||
  // When M is not specified, the identity is used.
  double maxDot = 0.0;

  int xc = (X) ? X->NumVectors() : 0;
  int rc = (R) ? R->NumVectors() : 0;
  
  if (xc*rc == 0)
    return maxDot;

  int i, j;
  for (i = 0; i < rc; ++i) {
    Epetra_Vector MRi(Epetra_DataAccess::Copy, (*R), i);
    if (M)
      M->Apply(*((*R)(i)), MRi);
    double normMR = 0.0;
    MRi.Norm2(&normMR);
    double dot = 0.0;
    for (j = 0; j < xc; ++j) {
      double normXj = 0.0;
      (*X)(j)->Norm2(&normXj);
      (*X)(j)->Dot(MRi, &dot);
      dot = fabs(dot)/(normMR*normXj);
      maxDot = (dot > maxDot) ? dot : maxDot;
    }
  }

  return maxDot;

}


double CheckingTools::errorOrthonormality(const Epetra_MultiVector *X, 
                      const Epetra_Operator *M) const {

  // Return the maximum coefficient of the matrix X^T * M * X - I
  // When M is not specified, the identity is used.
  double maxDot = 0.0;

  int xc = (X) ? X->NumVectors() : 0;
  if (xc == 0)
    return maxDot;

  int i, j;
  for (i = 0; i < xc; ++i) {
    Epetra_Vector MXi(Epetra_DataAccess::Copy, (*X), i);
    if (M)
      M->Apply(*((*X)(i)), MXi);
    double dot = 0.0;
    for (j = 0; j < xc; ++j) {
      (*X)(j)->Dot(MXi, &dot);
      dot = (i == j) ? fabs(dot - 1.0) : fabs(dot);
      maxDot = (dot > maxDot) ? dot : maxDot;
    }
  }

  return maxDot;

}


double CheckingTools::errorEquality(const Epetra_MultiVector *X,
                      const Epetra_MultiVector *MX, const Epetra_Operator *M) const {

  // Return the maximum coefficient of the matrix M * X - MX
  // scaled by the maximum coefficient of MX.
  // When M is not specified, the identity is used.

  double maxDiff = 0.0;

  int xc = (X) ? X->NumVectors() : 0;
  int mxc = (MX) ? MX->NumVectors() : 0;

  if ((xc != mxc) || (xc*mxc == 0))
    return maxDiff;

  int i;
  double maxCoeffX = 0.0;
  for (i = 0; i < xc; ++i) {
    double tmp = 0.0;
    (*MX)(i)->NormInf(&tmp);
    maxCoeffX = (tmp > maxCoeffX) ? tmp : maxCoeffX;
  }

  for (i = 0; i < xc; ++i) {
    Epetra_Vector MtimesXi(Epetra_DataAccess::Copy, (*X), i);
    if (M)
      M->Apply(*((*X)(i)), MtimesXi);
    MtimesXi.Update(-1.0, *((*MX)(i)), 1.0);
    double diff = 0.0;
    MtimesXi.NormInf(&diff);
    maxDiff = (diff > maxDiff) ? diff : maxDiff;
  }

  return (maxCoeffX == 0.0) ? maxDiff : maxDiff/maxCoeffX;

}


int CheckingTools::errorSubspaces(const Epetra_MultiVector &Q, const Epetra_MultiVector &Qex,
                                  const Epetra_Operator *M) const {

  using std::cout;
  using std::ios;

  int info = 0;
  int myPid = MyComm.MyPID();

  int qr = Q.MyLength();
  int qc = Q.NumVectors();
  int qexc = Qex.NumVectors();

  double *mQ = new (std::nothrow) double[qr];
  if (mQ == 0) {
    info = -1;
    return info;
  }

  double *z = new (std::nothrow) double[qexc*qc];
  if (z == 0) {
    delete[] mQ;
    info = -1;
    return info;
  }

  Epetra_LocalMap lMap(qexc, 0, MyComm);
  Epetra_MultiVector QextMQ(Epetra_DataAccess::View, lMap, z, qexc, qc);

  int j;
  for (j=0; j<qc; ++j) {
    Epetra_MultiVector Qj(Epetra_DataAccess::View, Q, j, 1);
    Epetra_MultiVector MQ(Epetra_DataAccess::View, Q.Map(), mQ, qr, 1);
    if (M)
      M->Apply(Qj, MQ);
    else
      memcpy(mQ, Qj.Values(), qr*sizeof(double));
    Epetra_MultiVector colJ(Epetra_DataAccess::View, QextMQ, j, 1);
    colJ.Multiply('T', 'N', 1.0, Qex, MQ,  0.0);
  }
  delete[] mQ;

  // Compute the SVD

  int svSize = (qc > qexc) ? qexc : qc;
  double *sv = new (std::nothrow) double[svSize];
  if (sv == 0) {
    delete[] z;
    info  = -1;
    return info;
  }

  // lwork is larger than the value suggested by LAPACK
  int lwork = (qexc > qc) ? qexc + 5*qc : qc + 5*qexc;
  double *work = new (std::nothrow) double[lwork];
  if (work == 0) {
    delete[] z;
    delete[] sv;
    info = -1;
    return info;
  }

  Epetra_LAPACK call;
  call.GESVD('N','N',qexc,qc,QextMQ.Values(),qexc,sv,0,qc,0,qc,work,&lwork,&info);

  delete[] work;
  delete[] z;

  // Treat error messages

  if (info < 0) {
    if (myPid == 0) {
      std::cerr << std::endl;
      std::cerr << " In DGESVD, argument " << -info << " has an illegal value\n";
      std::cerr << std::endl;
    }
    delete[] sv;
    return info;
  }

  if (info > 0) {
    if (myPid == 0) {
      std::cerr << std::endl;
      std::cerr << " In DGESVD, DBSQR did not converge (" << info << ").\n";
      std::cerr << std::endl;
    }
    delete[] sv;
    return info;
  }

  // Output the result

  if (myPid == 0) {
    std::cout << std::endl;
    std::cout << " --- Principal angles between eigenspaces ---\n";
    std::cout << std::endl;
    std::cout << " Reference space built with " << Qex.NumVectors() << " vectors." << std::endl;
    std::cout << " Dimension of computed space: " << Q.NumVectors() << std::endl;
    std::cout << std::endl;
    std::cout.setf(ios::scientific, ios::floatfield);
    std::cout.precision(4);
    std::cout << " Smallest singular value = " << sv[0] << std::endl;
    std::cout << " Largest singular value  = " << sv[svSize-1] << std::endl;
    std::cout << std::endl;
    std::cout << " Smallest angle between subspaces (rad) = ";
    std::cout << ((sv[0] > 1.0) ? 0.0 : asin(sqrt(1.0 - sv[0]*sv[0]))) << std::endl;
    std::cout << " Largest angle between subspaces (rad)  = ";
    std::cout << ((sv[0] > 1.0) ? 0.0 : asin(sqrt(1.0 - sv[svSize-1]*sv[svSize-1]))) << std::endl;
    std::cout << std::endl;
  }

  delete[] sv;

  return info;

}


void CheckingTools::errorEigenResiduals(const Epetra_MultiVector &Q, double *lambda,
                                  const Epetra_Operator *K, const Epetra_Operator *M,
                                  double *normWeight) const {

  using std::cout;
  using std::ios;

  if ((K == 0) || (lambda == 0))
    return;

  int myPid = MyComm.MyPID();

  int qr = Q.MyLength();
  int qc = Q.NumVectors();

  double *work = new (std::nothrow) double[2*qr];
  if (work == 0)
    return;

  if (myPid == 0) {
    std::cout << std::endl;
    std::cout << " --- Norms of residuals for computed eigenmodes ---\n";
    std::cout << std::endl;
    std::cout << "        Eigenvalue";
    if (normWeight)
      std::cout << "     User Norm     Scaled User N.";
    std::cout << "     2-Norm     Scaled 2-Nor.\n";
  }

  Epetra_Vector KQ(Epetra_DataAccess::View, Q.Map(), work);
  Epetra_Vector MQ(Epetra_DataAccess::View, Q.Map(), work + qr);
  Epetra_Vector Qj(Epetra_DataAccess::View, Q.Map(), Q.Values());

  double maxUserNorm = 0.0;
  double minUserNorm = 1.0e+100;
  double maxL2Norm = 0.0;
  double minL2Norm = 1.0e+100;

  // Compute the residuals and norms
  int j;
  for (j=0; j<qc; ++j) {

    Qj.ResetView(Q.Values() + j*qr);
    if (M)
      M->Apply(Qj, MQ);
    else
      memcpy(MQ.Values(), Qj.Values(), qr*sizeof(double));
    K->Apply(Qj, KQ);
    KQ.Update(-lambda[j], MQ, 1.0);

    double residualL2 = 0.0;
    KQ.Norm2(&residualL2);

    double residualUser = 0.0;
    if (normWeight) {
      Epetra_Vector vectWeight(Epetra_DataAccess::View, Q.Map(), normWeight);
      KQ.NormWeighted(vectWeight, &residualUser);
    }

    if (myPid == 0) {
      std::cout << " ";
      std::cout.width(4);
      std::cout.precision(8);
      std::cout.setf(ios::scientific, ios::floatfield);
      std::cout << j+1 << ". " << lambda[j] << " ";
      if (normWeight) {
        std::cout << residualUser << " ";
        std::cout << ((lambda[j] == 0.0) ? 0.0 : residualUser/lambda[j]) << " ";
      }
      std::cout << residualL2 << " " << ((lambda[j] == 0.0) ? 0.0 : residualL2/lambda[j]) << " ";
      std::cout << std::endl;
    }

    if (lambda[j] > 0.0) {
      maxL2Norm = (residualL2/lambda[j] > maxL2Norm) ? residualL2/lambda[j] : maxL2Norm;
      minL2Norm = (residualL2/lambda[j] < minL2Norm) ? residualL2/lambda[j] : minL2Norm;
      if (normWeight) {
        maxUserNorm = (residualUser/lambda[j] > maxUserNorm) ? residualUser/lambda[j]
                                                             : maxUserNorm;
        minUserNorm = (residualUser/lambda[j] < minUserNorm) ? residualUser/lambda[j]
                                                             : minUserNorm;
      }
    }

  } // for j=0; j<qc; ++j) 

  if (myPid == 0) {
    std::cout << std::endl;
    if (normWeight) {
      std::cout << " >> Minimum scaled user-norm of residuals = " << minUserNorm << std::endl;
      std::cout << " >> Maximum scaled user-norm of residuals = " << maxUserNorm << std::endl;
      std::cout << std::endl;
    }
    std::cout << " >> Minimum scaled L2 norm of residuals   = " << minL2Norm << std::endl;
    std::cout << " >> Maximum scaled L2 norm of residuals   = " << maxL2Norm << std::endl;
    std::cout << std::endl;
  }

  if (work)
    delete[] work;

}


void CheckingTools::errorEigenResiduals(const Epetra_MultiVector &Q, double *lambda,
                                  const Epetra_Operator *K, const Epetra_Operator *M,
                                  const Epetra_Operator *Msolver) const {

  using std::cout;
  using std::ios;

  if ((K == 0) || (lambda == 0) || (Msolver == 0))
    return;

  int myPid = MyComm.MyPID();

  int qr = Q.MyLength();
  int qc = Q.NumVectors();

  double *work = new (std::nothrow) double[2*qr];
  if (work == 0)
    return;

  if (myPid == 0) {
    std::cout << std::endl;
    std::cout << " --- Norms of residuals for computed eigenmodes ---\n";
    std::cout << std::endl;
    std::cout << "        Eigenvalue";
    std::cout << "     M^{-1} N.     Sca. M^{-1} N.";
    std::cout << "     2-Norm     Scaled 2-Nor.\n";
  }

  Epetra_Vector KQ(Epetra_DataAccess::View, Q.Map(), work);
  Epetra_Vector MQ(Epetra_DataAccess::View, Q.Map(), work + qr);
  Epetra_Vector Qj(Epetra_DataAccess::View, Q.Map(), Q.Values());

  double maxMinvNorm = 0.0;
  double minMinvNorm = 1.0e+100;
  double maxL2Norm = 0.0;
  double minL2Norm = 1.0e+100;

  // Compute the residuals and norms
  int j;
  for (j=0; j<qc; ++j) {

    Qj.ResetView(Q.Values() + j*qr);
    if (M)
      M->Apply(Qj, MQ);
    else
      memcpy(MQ.Values(), Qj.Values(), qr*sizeof(double));
    K->Apply(Qj, KQ);
    KQ.Update(-lambda[j], MQ, 1.0);

    double residualL2 = 0.0;
    KQ.Norm2(&residualL2);

    double residualMinv = 0.0;
    Msolver->ApplyInverse(KQ, MQ); 
    KQ.Dot(MQ, &residualMinv);
    residualMinv = sqrt(fabs(residualMinv));

    if (myPid == 0) {
      std::cout << " ";
      std::cout.width(4);
      std::cout.precision(8);
      std::cout.setf(ios::scientific, ios::floatfield);
      std::cout << j+1 << ". " << lambda[j] << " ";
      std::cout << residualMinv << " ";
      std::cout << ((lambda[j] == 0.0) ? 0.0 : residualMinv/lambda[j]) << " ";
      std::cout << residualL2 << " " << ((lambda[j] == 0.0) ? 0.0 : residualL2/lambda[j]) << " ";
      std::cout << std::endl;
    }

    if (lambda[j] > 0.0) {
      maxL2Norm = (residualL2/lambda[j] > maxL2Norm) ? residualL2/lambda[j] : maxL2Norm;
      minL2Norm = (residualL2/lambda[j] < minL2Norm) ? residualL2/lambda[j] : minL2Norm;
      maxMinvNorm = (residualMinv/lambda[j] > maxMinvNorm) ? residualMinv/lambda[j]
                                                           : maxMinvNorm;
      minMinvNorm = (residualMinv/lambda[j] < minMinvNorm) ? residualMinv/lambda[j]
                                                           : minMinvNorm;
    }

  } // for j=0; j<qc; ++j) 

  if (myPid == 0) {
    std::cout << std::endl;
    std::cout << " >> Minimum scaled M^{-1}-norm of residuals = " << minMinvNorm << std::endl;
    std::cout << " >> Maximum scaled M^{-1}-norm of residuals = " << maxMinvNorm << std::endl;
    std::cout << std::endl;
    std::cout << " >> Minimum scaled L2 norm of residuals   = " << minL2Norm << std::endl;
    std::cout << " >> Maximum scaled L2 norm of residuals   = " << maxL2Norm << std::endl;
    std::cout << std::endl;
  }

  if (work)
    delete[] work;

}


int CheckingTools::errorLambda(double *continuous, double *discrete, int numDiscrete,
                               double *lambda, int nev, bool smallest) const {

  using std::cout;
  using std::ios;

  int myPid = MyComm.MyPID();
  int nMax = 0;

  // Allocate working arrays

  int *used = new (std::nothrow) int[numDiscrete + nev];
  if (used == 0) {
    return 0;
  }

  int *bestMatch = used + numDiscrete;

  // Find the best match for the eigenvalues
  double eps = 0.0;
  Epetra_LAPACK call;
  call.LAMCH('E', eps);

  double gap = Epetra_MaxDouble;
  for (int i=0; i<numDiscrete; ++i) {
    used[i] = -1;
    for (int j = i; j < numDiscrete; ++j)
      if (discrete[j] > (1.0 + 10.0*eps)*discrete[i]) {
        double tmp = (discrete[j] - discrete[i])/discrete[i];
        gap = (tmp < gap) ? tmp : gap;
        break;
      }
  }

  for (int i=0; i<nev; ++i) 
    bestMatch[i] = -1;

  for (int i=0; i<nev; ++i) {
    if (lambda[i] < continuous[0])
      continue;
    bestMatch[i] = (i == 0) ? 0 : bestMatch[i-1] + 1;
    int jStart = bestMatch[i];
    for (int j = jStart; j < numDiscrete; ++j) {
      double diff = fabs(lambda[i]-discrete[j]);
      if (diff < 0.5*gap*lambda[i]) {
        bestMatch[i] = j;
        break;
      }
    }
    used[bestMatch[i]] = i;
  }

  // Print the results for the eigenvalues
  if (myPid == 0) {
    std::cout << std::endl;
    std::cout << " --- Relative errors on eigenvalues ---\n";
    std::cout << std::endl;
    std::cout << "       Exact Cont.    Exact Disc.     Computed      ";
    std::cout << " Alg. Err.  Mesh Err.\n";
  }

  int iCount = 0;
  for (int i=0; i<nev; ++i) {
    if (bestMatch[i] == -1) {
      if (myPid == 0) {
        std::cout << "      ************** ************** ";
        std::cout.precision(8);
        std::cout.setf(ios::scientific, ios::floatfield);
        std::cout << lambda[i];
        std::cout << "  *********  *********" << std::endl;
      }
      iCount += 1;
    }
  }

  double lastDiscrete = 0.0;
  if (smallest) {
    // print the smallest eigenvalues, up to the ones we've matched
    for (int i=0; i<numDiscrete; ++i) {
      if ((iCount == nev) && (discrete[i] > lastDiscrete))
        break;
      nMax = i;
      if (used[i] < 0) {
        lastDiscrete = discrete[i];
        if (myPid == 0) {
          std::cout << " ";
          std::cout.width(4);
          std::cout << i+1 << ". ";
          std::cout.precision(8);
          std::cout.setf(ios::scientific, ios::floatfield);
          std::cout << continuous[i] << " " << discrete[i] << " ";
          std::cout << "**************  *********  ";
          std::cout.precision(3);
          std::cout << fabs(continuous[i]-discrete[i])/continuous[i] << std::endl;
        }
      }
      else {
        lastDiscrete = discrete[i];
        if (myPid == 0) {
          std::cout << " ";
          std::cout.width(4);
          std::cout << i+1 << ". ";
          std::cout.precision(8);
          std::cout.setf(ios::scientific, ios::floatfield);
          std::cout << continuous[i] << " " << discrete[i] << " " << lambda[used[i]] << "  ";
          std::cout.precision(3);
          std::cout << fabs(lambda[used[i]]-discrete[i])/discrete[i] << "  ";
          std::cout << fabs(continuous[i]-discrete[i])/continuous[i] << std::endl;
        }
        iCount += 1;
      }
    }
  }
  else {
    // print the largest eigenvalues, from the first that we've matched
    bool started = false;
    for (int i=0; i<numDiscrete; ++i) {
      if (used[i] >= 0) {
        if (started == false) {
          nMax = i;
        }
        started = true;
        if (myPid == 0) {
          std::cout << " ";
          std::cout.width(4);
          std::cout << i+1 << ". ";
          std::cout.precision(8);
          std::cout.setf(ios::scientific, ios::floatfield);
          std::cout << continuous[i] << " " << discrete[i] << " " << lambda[used[i]] << "  ";
          std::cout.precision(3);
          std::cout << fabs(lambda[used[i]]-discrete[i])/discrete[i] << "  ";
          std::cout << fabs(continuous[i]-discrete[i])/continuous[i] << std::endl;
        }
        iCount += 1;
      }
      else if (started) {
        if (myPid == 0) {
          std::cout << " ";
          std::cout.width(4);
          std::cout << i+1 << ". ";
          std::cout.precision(8);
          std::cout.setf(ios::scientific, ios::floatfield);
          std::cout << continuous[i] << " " << discrete[i] << " ";
          std::cout << "**************  *********  ";
          std::cout.precision(3);
          std::cout << fabs(continuous[i]-discrete[i])/continuous[i] << std::endl;
        }
      }
    }
  }

  delete[] used;

  return nMax;

}


int CheckingTools::inputArguments(const int &numEigen, const Epetra_Operator *K, 
                            const Epetra_Operator *M, const Epetra_Operator *P,
                            const Epetra_MultiVector &Q, const int &minSize) const {

  using std::cout;
  using std::ios;

  // Routine to check some arguments
  // 
  // info = -  1 >> The stiffness matrix K has not been specified.
  // info = -  2 >> The maps for the matrix K and the matrix M differ.
  // info = -  3 >> The maps for the matrix K and the preconditioner P differ.
  // info = -  4 >> The maps for the vectors and the matrix K differ.
  // info = -  5 >> Q is too small for the number of eigenvalues requested.
  // info = -  6 >> Q is too small for the computation parameters.
  //

  int myPid = MyComm.MyPID();

  if (K == 0) {
    if (myPid == 0)
      std::cerr << "\n !!! The matrix K to analyze has not been specified !!! \n\n";
    return -1;
  }

  if (M) {
    int mGlobal = (M->OperatorDomainMap()).NumGlobalElements();
    int mLocal = (M->OperatorDomainMap()).NumMyElements();
    int kGlobal = (K->OperatorDomainMap()).NumGlobalElements();
    int kLocal = (K->OperatorDomainMap()).NumMyElements();
    if ((mGlobal != kGlobal) || (mLocal != kLocal)) {
      if (myPid == 0) {
        std::cerr << std::endl;
        std::cerr << " !!! The maps for the matrice K and the mass M are different !!!\n";
        std::cerr << std::endl;
      }
      return -2;
    }
  }

  if (P) {
    int pGlobal = (P->OperatorDomainMap()).NumGlobalElements();
    int pLocal = (P->OperatorDomainMap()).NumMyElements();
    int kGlobal = (K->OperatorDomainMap()).NumGlobalElements();
    int kLocal = (K->OperatorDomainMap()).NumMyElements();
    if ((pGlobal != kGlobal) || (pLocal != kLocal)) {
      if (myPid == 0) {
        std::cerr << std::endl;
        std::cerr << " !!! The maps for the matrice K and the preconditioner P are different !!!\n";
        std::cerr << std::endl;
      }
      return -3;
    }
  }

  if ((Q.MyLength() != (K->OperatorDomainMap()).NumMyElements()) || 
      (Q.GlobalLength() != (K->OperatorDomainMap()).NumGlobalElements())) {
    if (myPid == 0) {
      std::cerr << "\n !!! The maps for the vectors and the matrix are different !!! \n\n";
    }
    return -4;
  }

  if (Q.NumVectors() < numEigen) {
    if (myPid == 0) {
      std::cerr << std::endl;
      std::cerr << " !!! The number of eigenvalues is too large for the space allocated !!! \n\n";
      std::cerr << " The recommended size for " << numEigen << " eigenvalues is ";
      std::cerr << minSize << std::endl;
      std::cerr << std::endl;
    }
    return -5;
  }

  if (Q.NumVectors() < minSize) {
    if (myPid == 0) {
      std::cerr << std::endl;
      std::cerr << " !!! The space allocated is too small for the number of eigenvalues";
      std::cerr << " and the size of blocks specified !!! \n";
      std::cerr << " The recommended size for " << numEigen << " eigenvalues is ";
      std::cerr << minSize << std::endl;
      std::cerr << std::endl;
    }
    return -6;
  }

  return 0;

}
