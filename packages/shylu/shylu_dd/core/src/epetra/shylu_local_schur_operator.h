// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SHYLU_LOCAL_SCHUR_OPERATOR_H
#define SHYLU_LOCAL_SCHUR_OPERATOR_H

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Amesos_BaseSolver.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Teuchos_Time.hpp>
#include <Teuchos_RCP.hpp>
#include "shylu_util.h"
#include "shylu.h"
#include <EpetraExt_Reindex_MultiVector.h>

class ShyLU_Local_Schur_Operator : public virtual Epetra_Operator
{

    public:

    // TODO: Change to RCPs
    ShyLU_Local_Schur_Operator(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver,
    Ifpack_Preconditioner *ifSolver, Epetra_CrsMatrix *C,
    Epetra_Map *LocalDRowMap, int nvectors);

    ~ShyLU_Local_Schur_Operator()
    {
        ssym_->OrigLP->SetLHS(orig_lhs_);
        ssym_->OrigLP->SetRHS(orig_rhs_);
        ssym_->ReIdx_LP->fwd();
    }

    int SetUseTranspose(bool useTranspose);

    int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    double NormInf() const;

    const char *Label() const;

    bool UseTranspose() const;

    bool HasNormInf() const;

    const Epetra_Comm& Comm() const;

    const Epetra_Map& OperatorDomainMap() const;

    const Epetra_Map& OperatorRangeMap() const;

    void PrintTimingInfo();

    void ResetTempVectors(int nvectors);

    mutable int cntApply;

    Epetra_CrsMatrix *G_;
    Epetra_CrsMatrix *R_;
    Epetra_LinearProblem *LP_;
    Amesos_BaseSolver *solver_;
    Ifpack_Preconditioner *ifSolver_;
    Epetra_CrsMatrix *C_;
    Epetra_MultiVector *orig_lhs_;
    Epetra_MultiVector *orig_rhs_;
    Epetra_Map *localDRowMap_;

    int nvectors_;
    Teuchos::RCP<Epetra_MultiVector> temp;
    Teuchos::RCP<Epetra_MultiVector> temp2;
    Teuchos::RCP<Epetra_MultiVector> localX;
    shylu_symbolic *ssym_;   // symbolic structure
    shylu_config *config_;

#ifdef TIMING_OUTPUT
    Teuchos::RCP<Teuchos::Time> matvec_time_;
    Teuchos::RCP<Teuchos::Time> localize_time_;
    Teuchos::RCP<Teuchos::Time> trisolve_time_;
    Teuchos::RCP<Teuchos::Time> dist_time_;
    Teuchos::RCP<Teuchos::Time> matvec2_time_;
    Teuchos::RCP<Teuchos::Time> apply_time_;
    Teuchos::RCP<Teuchos::Time> update_time_;
#endif

};
#endif // SHYLU_LOCAL_SCHUR_OPERATOR_H
