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

#ifndef AMG_OPERATOR_H
#define AMG_OPERATOR_H

#include "Epetra_ConfigDefs.h"
#include "AztecOO.h"

#include "ml_include.h"
#include "ml_epetra_operator.h"
#include "ml_epetra_utils.h"
#include "ml_agg_METIS.h"

#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#ifndef MACOSX
#include "singularCoarse.h"
#endif

class AMGOperator : public virtual Epetra_Operator {

  private:

    const Epetra_Comm &MyComm;
    const Epetra_BLAS callBLAS;
    const Epetra_LAPACK callLAPACK;

    const Epetra_Operator *K;
    Epetra_Operator *Prec;

    const Epetra_MultiVector *Q;
    double *QtQ;

    int numDofs;

    bool leftProjection;
    bool rightProjection;

    ML *ml_handle;
    ML_Aggregate *ml_agg;

    int AMG_NLevels;

    int coarseLocalSize;
    int coarseGlobalSize;

    double *ZcoarseTZcoarse;

    int verbose;

    void preProcess(int maxCoarseSize);
    void setCoarseSolver_Cycle(int coarseSolver, int cycle);

    // Don't define these functions
    AMGOperator(const AMGOperator &ref);
    AMGOperator& operator=(const AMGOperator &ref);

  public:

    AMGOperator(const Epetra_Comm& _Com, const Epetra_Operator *KK, int verb = 0,
                int nLevel = 10, int smoother = 1, int param = 2,
                int coarseSolver = -1, int cycle = 0,
                int _numDofs = 1, const Epetra_MultiVector *Z = 0);

    AMGOperator(const Epetra_Comm& _Com, const Epetra_Operator *KK, int verb = 0,
                int nLevel = 10, int smoother = 1, int *param = 0,
                int coarseSolver = -1, int cycle = 0,
                int _numDofs = 1, const Epetra_MultiVector *Z = 0);

    int SetUseLeftProjection(bool proj) { leftProjection = proj; return 0; }
    int SetUseRightProjection(bool proj) { rightProjection = proj; return 0; }

    ~AMGOperator();

    char * Label() const { return "Epetra_Operator for AMG preconditioner"; };

    bool UseTranspose() const { return (false); };
    int SetUseTranspose(bool UseTranspose) { return 0; };

    bool HasNormInf() const { return (false); };
    double NormInf() const  { return (-1.0); };

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    const Epetra_Comm& Comm() const { return MyComm; };

    const Epetra_Map& OperatorDomainMap() const { return K->OperatorDomainMap(); };
    const Epetra_Map& OperatorRangeMap() const { return K->OperatorRangeMap(); };

    int getAMG_NLevels() const { return AMG_NLevels; };

};

#endif
