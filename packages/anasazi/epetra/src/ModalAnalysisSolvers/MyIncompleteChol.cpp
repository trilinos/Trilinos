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

#include "MyIncompleteChol.h"


MyIncompleteChol::MyIncompleteChol(const Epetra_Comm &_Comm, const Epetra_Operator *KK,
                         double tolDrop, int fillLevel, const Epetra_MultiVector *Z)
               : MyComm(_Comm),
                 callBLAS(),
                 callLAPACK(),
                 K(KK),
                 Prec(0),
                 dropTol(tolDrop),
                 lFill(fillLevel),
                 Q(Z),
                 QtQ(0),
                 leftProjection(false),
                 rightProjection(false)
               {

  // Note the constness is cast away for ML
  Epetra_CrsMatrix *Kcrs = dynamic_cast<Epetra_CrsMatrix*>(const_cast<Epetra_Operator*>(K));
  if (Kcrs == 0) {
    if (MyComm.MyPID() == 0) {
      cerr << endl;
      cerr << " !!! For incomplete Cholesky preconditioner, ";
      cerr << "the matrix must be 'Epetra_CrsMatrix' object !!!\n";
      cerr << endl;
    }
    return;
  }

  Prec = new Ifpack_CrsIct(*Kcrs, dropTol, Kcrs->GlobalMaxNumEntries() + lFill);
  Prec->InitValues(*Kcrs);
  Prec->Factor();

  if (Q) {
    QtQ = new double[Q->NumVectors()*Q->NumVectors()];
    double *work = new double[Q->NumVectors()*Q->NumVectors()];
    callBLAS.GEMM('T', 'N', Q->NumVectors(), Q->NumVectors(), Q->MyLength(),
                  1.0, Q->Values(), Q->MyLength(), Q->Values(), Q->MyLength(),
                  0.0, work, Q->NumVectors());
    MyComm.SumAll(work, QtQ, Q->NumVectors()*Q->NumVectors());
    delete[] work;
    int info = 0;
    callLAPACK.POTRF('U', Q->NumVectors(), QtQ, Q->NumVectors(), &info);
    if (info) {
      cerr << endl;
      cerr << " !!! In incomplete Cholesky, the null space vectors are linearly dependent !!!\n";
      cerr << endl;
      delete[] QtQ;
      QtQ = 0;
    }
  }

}


MyIncompleteChol::~MyIncompleteChol() {

  if (Prec) {
    delete Prec;
    Prec = 0;
  }

  if (QtQ) {
    delete[] QtQ;
    QtQ = 0;
  }

}


int MyIncompleteChol::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  Y.PutScalar(0.0);

  return -1;

}


int MyIncompleteChol::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  int xcol = X.NumVectors();

  if (Y.NumVectors() < xcol)
    return -1;

  double *valX = (rightProjection == true) ? new double[xcol*X.MyLength()] : X.Values();

  Epetra_MultiVector X2(View, X.Map(), valX, X.MyLength(), xcol);

  if ((rightProjection == true) && (Q)) {
    int qcol = Q->NumVectors();
    double *work = new double[2*qcol*xcol];
    memcpy(X2.Values(), X.Values(), xcol*X.MyLength()*sizeof(double)); 
    callBLAS.GEMM('T', 'N', qcol, xcol, Q->MyLength(), 1.0, Q->Values(), Q->MyLength(),
                  X2.Values(), X2.MyLength(), 0.0, work + qcol*xcol, qcol);
    MyComm.SumAll(work + qcol*xcol, work, qcol*xcol);
    int info = 0;
    callLAPACK.POTRS('U', qcol, xcol, QtQ, qcol, work, qcol, &info); 
    callBLAS.GEMM('N', 'N', X2.MyLength(), xcol, qcol, -1.0, Q->Values(), Q->MyLength(),
                  work, qcol, 1.0, X2.Values(), X2.MyLength());
    delete[] work;
  }

  Prec->ApplyInverse(X2, Y);

  if (rightProjection)
    delete[] valX;

  if ((leftProjection == true) && (Q)) {
    int qcol = Q->NumVectors();
    double *work = new double[2*qcol*xcol];
    callBLAS.GEMM('T', 'N', qcol, xcol, Q->MyLength(), 1.0, Q->Values(), Q->MyLength(),
                  Y.Values(), Y.MyLength(), 0.0, work + qcol*xcol, qcol);
    MyComm.SumAll(work + qcol*xcol, work, qcol*xcol);
    int info = 0;
    callLAPACK.POTRS('U', qcol, xcol, QtQ, qcol, work, qcol, &info); 
    callBLAS.GEMM('N', 'N', Y.MyLength(), xcol, qcol, -1.0, Q->Values(), Q->MyLength(),
                  work, qcol, 1.0, Y.Values(), Y.MyLength());
    delete[] work;
  }

  return 0;

}


