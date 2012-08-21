
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


#ifndef SHYLU_LOCAL_SCHUR_OPERATOR_H
#define SHYLU_LOCAL_SCHUR_OPERATOR_H

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
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
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
    Epetra_CrsMatrix *C_;
    Epetra_MultiVector *orig_lhs_;
    Epetra_MultiVector *orig_rhs_;
    Epetra_Map *localDRowMap_;

    int nvectors_;
    Teuchos::RCP<Epetra_MultiVector> temp;
    Teuchos::RCP<Epetra_MultiVector> temp2;
    Teuchos::RCP<Epetra_MultiVector> localX;
    shylu_symbolic *ssym_;   // symbolic structure

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
