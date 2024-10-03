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

#ifndef MY_INCOMPLETE_CHOL_H
#define MY_INCOMPLETE_CHOL_H

#include "Epetra_ConfigDefs.h"
#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_CrsIct.h"   

class MyIncompleteChol : public virtual Epetra_Operator {

  private:

    const Epetra_Comm &MyComm;
    const Epetra_BLAS callBLAS;
    const Epetra_LAPACK callLAPACK;

    const Epetra_Operator *K;

    Ifpack_CrsIct *Prec;
    double dropTol;
    int lFill;

    const Epetra_MultiVector *Q;
    double *QtQ;

    bool leftProjection;
    bool rightProjection;

    // Don't define these functions
    MyIncompleteChol(const MyIncompleteChol &ref);
    MyIncompleteChol& operator=(const MyIncompleteChol &ref);

  public:

    MyIncompleteChol(const Epetra_Comm& _Com, const Epetra_Operator *KK,
                double dropTol = Epetra_MinDouble, int lFill = 0, const Epetra_MultiVector *Z = 0);

    int SetUseLeftProjection(bool proj) { leftProjection = proj; return 0; }
    int SetUseRightProjection(bool proj) { rightProjection = proj; return 0; }

    ~MyIncompleteChol();

    char * Label() const { return "Epetra_Operator for incomplete Cholesky preconditioner"; };

    bool UseTranspose() const { return (false); };
    int SetUseTranspose(bool UseTranspose) { return 0; };

    bool HasNormInf() const { return (false); };
    double NormInf() const  { return (-1.0); };

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    const Epetra_Comm& Comm() const { return MyComm; };

    const Epetra_Map& OperatorDomainMap() const { return K->OperatorDomainMap(); };
    const Epetra_Map& OperatorRangeMap() const { return K->OperatorRangeMap(); };

};

#endif
